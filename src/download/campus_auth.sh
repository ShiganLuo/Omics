#!/usr/bin/env bash
# campus_auth.sh — Nankai campus network (Dr.COM eportal) keep-alive + auto login.
#
# The campus gateway intercepts all HTTP/HTTPS (fake "Web secure CA" cert,
# DNS forgery) until the client is portal-authenticated. This script checks
# reachability and re-logs-in via the Dr.COM eportal API when needed.
#
# Usage:
#   campus_auth.sh [OPTIONS] [COMMAND]
#
# Commands (default: check + login if offline):
#   --status            only report status, never login
#   --login             force login even if probe says online
#   --logout            portal logout
#
# Options:
#   -u, --user USER     campus account (overrides ~/.campus_auth)
#   -p, --pass PASS     campus password (overrides ~/.campus_auth;
#                       if omitted at login time, prompted interactively)
#   -h, --help          show this help
#
# Credential precedence: -u/-p > ~/.campus_auth (chmod 600) > interactive
# prompt. When -u/-p are given, they are saved to ~/.campus_auth so the
# cron keep-alive can reuse them.
#
# Examples:
#   campus_auth.sh --login -u 1234567 -p mypass   # login now + persist creds
#   campus_auth.sh                                # cron mode: fix only if offline
#   campus_auth.sh --status
#
# Cron (every 10 min):
#   */10 * * * * /data/pub/zhousha/Totipotent20251031/workflow/Omics/src/download/campus_auth.sh >> /home/zhousha/.campus_auth.log 2>&1

set -u

# cron runs with a minimal env; ensure HOME resolves even there
: "${HOME:=$(getent passwd "$(id -un)" | cut -d: -f6)}"

EPORTAL="http://netauth.nankai.edu.cn:801/eportal"
PROBE_URL="https://www.baidu.com"
CRED_FILE="${HOME}/.campus_auth"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

usage() { sed -n '2,35p' "$0" | sed 's/^# \{0,1\}//'; }

my_ip() {
    # IP as seen by the campus gateway (source IP on the default route)
    ip -4 route get 223.5.5.5 2>/dev/null | awk '{for (i=1;i<=NF;i++) if ($i=="src") {print $(i+1); exit}}'
}

online() {
    local code
    code=$(curl -s --noproxy '*' --max-time 8 -o /dev/null -w '%{http_code}' "$PROBE_URL" 2>/dev/null)
    [[ "$code" == "200" ]]
}

load_creds_from_file() {
    [[ -r "$CRED_FILE" ]] || return 1
    # shellcheck source=/dev/null
    source "$CRED_FILE"
    [[ -n "${USER:-}" && -n "${PASS:-}" ]]
}

save_creds_to_file() {
    ( umask 077
      printf 'USER=%q\nPASS=%q\n' "$1" "$2" > "$CRED_FILE" )
    log "Credentials saved to $CRED_FILE (chmod 600)"
}

resolve_creds() {
    # fills CRED_USER / CRED_PASS: CLI args > cred file > interactive prompt
    CRED_USER="${ARG_USER:-}"; CRED_PASS="${ARG_PASS:-}"
    if [[ -n "$CRED_USER" && -n "$CRED_PASS" ]]; then
        save_creds_to_file "$CRED_USER" "$CRED_PASS"
        return 0
    fi
    if load_creds_from_file; then
        CRED_USER="$USER"; CRED_PASS="$PASS"
        return 0
    fi
    if [[ -z "$CRED_USER" ]]; then
        if [[ -t 0 ]]; then
            printf 'Campus account: '; read -r CRED_USER
        else
            log "ERROR: no account. Use -u/--user, or create $CRED_FILE (USER=/PASS= lines)."
            return 1
        fi
    fi
    if [[ -z "$CRED_PASS" ]]; then
        if [[ -t 0 ]]; then
            printf 'Campus password: '; read -rs CRED_PASS; printf '\n'
        else
            log "ERROR: no password. Use -p/--pass, or create $CRED_FILE."
            return 1
        fi
    fi
    save_creds_to_file "$CRED_USER" "$CRED_PASS"
}

do_login() {
    resolve_creds || return 1

    local ip resp
    ip=$(my_ip)
    if [[ -z "$ip" ]]; then
        log "ERROR: cannot determine source IP"
        return 1
    fi

    resp=$(curl -sk --noproxy '*' --max-time 10 \
        "${EPORTAL}/?c=ACSetting&a=Login&url=drappall" \
        --data-urlencode "DDDDD=,0,${CRED_USER}" \
        --data-urlencode "upass=${CRED_PASS}" \
        --data "0MKKey=123456789&v6ip=&wlanacname=blt2_&wlanuserip=${ip}" 2>/dev/null)

    if grep -q 'Dr.COMWebLoginID_3.htm' <<<"$resp"; then
        log "Login succeed (account=${CRED_USER}, ip=${ip})"
        return 0
    elif grep -q 'Dr.COMWebLoginID_2.htm' <<<"$resp"; then
        log "Login rejected by portal (bad credentials, or already online — check with --status)"
        return 1
    else
        log "Login FAILED: unexpected response: ${resp:0:200}"
        return 1
    fi
}

do_logout() {
    curl -sk --noproxy '*' --max-time 10 \
        "${EPORTAL}/?c=ACSetting&a=Logout&url=drappall" -o /dev/null
    log "Logout sent"
}

# ---- argument parsing (options and command may appear in any order) ----
COMMAND="check"
ARG_USER=""; ARG_PASS=""
while (($#)); do
    case "$1" in
        --status)  COMMAND="status" ;;
        --login)   COMMAND="login" ;;
        --logout)  COMMAND="logout" ;;
        -u|--user) ARG_USER="${2:?--user needs a value}"; shift ;;
        -p|--pass) ARG_PASS="${2:?--pass needs a value}"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log "Unknown argument: $1"; usage; exit 2 ;;
    esac
    shift
done

case "$COMMAND" in
    status)
        if online; then log "ONLINE"; else log "OFFLINE (gateway intercepting)"; fi
        ;;
    logout)
        do_logout
        ;;
    login)
        do_login
        ;;
    check)
        if online; then
            log "ONLINE"
        else
            log "OFFLINE, attempting portal login..."
            do_login
        fi
        ;;
esac
