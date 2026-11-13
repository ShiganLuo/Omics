#!/usr/bin/env python3
"""Execution-completion callback for the BI platform.

Posts one JSON message to the platform's public execution callback endpoint
when a workflow run finishes (SUCCESS), dies (FAILED) or is terminated
(CANCELLED). The endpoint is login-free and authenticated per execution by
a dedicated token:

    POST http://127.0.0.1:8083/api/public/executions/callback
    {
      "executionId": 12,
      "token": "<per-execution token>",
      "status": "SUCCESS | FAILED | CANCELLED",
      "message": "failure reason etc. (optional)",
      "exitCode": 0
    }

Design rules: exactly one send per ExecutionCallback instance (double-send
protection for signal-vs-exception races), never raise — a dead callback
endpoint must not change the workflow's exit code — bounded timeout with
one retry.

Usage:
    cb = ExecutionCallback(url, execution_id, token)
    cb.send("SUCCESS", exit_code=0)
    cb.send("FAILED", message="RuntimeError: ...", exit_code=1)
    cb.send("CANCELLED", message="terminated by signal 15", exit_code=143)
"""

import json
import logging
import urllib.error
import urllib.request
from typing import Any, Dict, Optional

DEFAULT_CALLBACK_URL = "http://127.0.0.1:8083/api/public/executions/callback"


class ExecutionCallback:
    """Once-guarded HTTP POST of an execution status to the callback endpoint.

    Parameters
    ----------
    url : str
        Callback endpoint (POST, JSON body).
    execution_id : str or int
        Platform execution id sent as ``executionId``; numeric strings are
        coerced to int to match the API contract.
    token : str
        Per-execution token sent as ``token``.
    timeout : float, default 10
        Socket timeout in seconds per attempt.
    attempts : int, default 2
        Total send attempts (first plus one retry).
    logger : logging.Logger, optional
        Logger for send diagnostics; a module-level logger is used when
        omitted.

    Attributes
    ----------
    sent : bool
        True after the first send() call; later calls are no-ops.
    """

    def __init__(self, url: str, execution_id: Any, token: str,
                 timeout: float = 10.0, attempts: int = 2,
                 logger: Optional[logging.Logger] = None) -> None:
        self.url = url
        try:
            self.execution_id: Any = int(execution_id)
        except (TypeError, ValueError):
            self.execution_id = execution_id
        self.token = token
        self.timeout = timeout
        self.attempts = max(1, attempts)
        self.logger = logger or logging.getLogger("execution_callback")
        self.sent = False

    def build_payload(self, status: str, message: str = "",
                      exit_code: int = 0) -> Dict[str, Any]:
        """Build the JSON-serializable callback body.

        Parameters
        ----------
        status : str
            One of SUCCESS / FAILED / CANCELLED.
        message : str, optional
            Human-readable detail (failure reason etc.); omitted when empty.
        exit_code : int, optional
            Process exit code, recorded only.

        Returns
        -------
        dict
            Payload with executionId, token, status and the optional
            message / exitCode fields.
        """
        payload: Dict[str, Any] = {
            "executionId": self.execution_id,
            "token": self.token,
            "status": status,
        }
        if message:
            payload["message"] = message
        if exit_code is not None:
            payload["exitCode"] = int(exit_code)
        return payload

    def send(self, status: str, message: str = "", exit_code: int = 0) -> bool:
        """POST the callback once; later calls are ignored.

        Returns
        -------
        bool
            True when the endpoint answered 2xx; False on any failure
            (logged as a warning — the caller's exit path is unaffected).
        """
        if self.sent:
            return False
        self.sent = True
        payload = self.build_payload(status, message, exit_code)
        body = json.dumps(payload).encode("utf-8")
        last_err = ""
        for attempt in range(1, self.attempts + 1):
            req = urllib.request.Request(
                self.url, data=body, method="POST",
                headers={"Content-Type": "application/json"})
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    code = getattr(resp, "status", None) or resp.getcode()
                    if 200 <= int(code) < 300:
                        self.logger.info(
                            f"execution callback sent: status={status} "
                            f"exitCode={exit_code} http={code}")
                        return True
                    last_err = f"HTTP {code}"
            except (urllib.error.URLError, urllib.error.HTTPError,
                    OSError, ValueError) as e:
                last_err = f"{type(e).__name__}: {e}"
            self.logger.warning(
                f"execution callback attempt {attempt}/{self.attempts} "
                f"failed ({last_err})")
        self.logger.warning(
            f"execution callback dropped after {self.attempts} attempts: "
            f"{self.url}")
        return False
