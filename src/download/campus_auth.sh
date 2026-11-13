#!/usr/bin/env bash
# campus_auth.sh — campus gateway keep-alive + auto login (multi-gateway).
#
# The campus gateway intercepts all outbound HTTP/HTTPS (302 redirect to a
# portal, or fake cert on TLS) until the client is portal-authenticated.
# This script probes reachability and re-logs-in via the matching portal
# API when needed. Designed to run from cron for unattended keep-alive.
#
# Gateways, auto-detected from the intercept redirect:
#   panabit : Panabit webauth portal (lv_zheda network; oauth2.*:8080)
#             POST /api?route=webauth&action=user_login
#             password AES-128-ECB/ZeroPadding/hex (key: Panabit@1024_key)
#   drcom   : Dr.COM eportal (Nankai; netauth.nankai.edu.cn:801)
#             POST /eportal?c=ACSetting&a=Login
#
# Usage:
#   campus_auth.sh [OPTIONS] [COMMAND]
#
# Commands (default: check + login if offline):
#   --status            only report status, never login
#   --login             force login even if probe says online
#   --logout            portal logout (drcom only; panabit has no logout API)
#   --selftest          verify the AES implementation against known vectors
#
# Options:
#   -u, --user USER     campus account (overrides ~/.campus_auth)
#   -p, --pass PASS     campus password (overrides ~/.campus_auth;
#                       if omitted at login time, prompted interactively)
#   -h, --help          show this help
#
# Env:
#   PROBE_URL     probe URL (default https://www.baidu.com); point at a
#                 dead host to simulate "offline" and exercise the login path
#   PA_PORTAL     Panabit portal base (default http://oauth2.panabit.com:8080)
#
# Credential precedence: -u/-p > ~/.campus_auth (chmod 600) > interactive
# prompt. When -u/-p are given, they are saved to ~/.campus_auth so the
# cron keep-alive can reuse them. After 3 consecutive failures the login
# is paused 30 min to avoid tripping the portal's account lockout.
#
# Examples:
#   campus_auth.sh --login -u 1234567 -p mypass   # login now + persist creds
#   campus_auth.sh                                # cron mode: fix only if offline
#   campus_auth.sh --status
#
# Cron (every 5 min; silent while healthy):
#   */5 * * * * /data/16T/luosg/workflow/Omics/src/download/campus_auth.sh >> /home/luosg/.campus_auth.log 2>&1

set -u

# cron runs with a minimal env; ensure HOME resolves even there
: "${HOME:=$(getent passwd "$(id -un)" | cut -d: -f6)}"

PROBE_URL="${PROBE_URL:-https://www.baidu.com}"
PA_PORTAL="${PA_PORTAL:-http://oauth2.panabit.com:8080}"
EPORTAL="http://netauth.nankai.edu.cn:801/eportal"
CRED_FILE="${HOME}/.campus_auth"
STATE_FILE="${HOME}/.campus_auth.state"
PA_KEY_HEX="50616e6162697440313032345f6b6579"   # AES-128 key "Panabit@1024_key"
FAIL_PAUSE=1800                                 # pause after 3 login failures

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

usage() { awk 'NR>1 && !/^#/ {exit} NR>1 {sub(/^# ?/,""); print}' "$0"; }

my_ip() {
    # IP as seen by the campus gateway (source IP on the default route)
    ip -4 route get 223.5.5.5 2>/dev/null | awk '{for (i=1;i<=NF;i++) if ($i=="src") {print $(i+1); exit}}'
}

my_mac() {
    local dev
    dev=$(ip -4 route get 223.5.5.5 2>/dev/null | grep -oP 'dev \K\S+' | head -1)
    [[ -n "$dev" ]] && cat "/sys/class/net/${dev}/address" 2>/dev/null
}

online() {
    local code url
    for url in "$PROBE_URL" https://www.qq.com; do
        code=$(curl -s --noproxy '*' --max-time 8 -o /dev/null -w '%{http_code}' "$url" 2>/dev/null)
        [[ "$code" == "200" ]] && return 0
    done
    return 1
}

# ---------- credentials ----------
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

# failure cooldown: skip login attempts for FAIL_PAUSE after 3 consecutive
# failures, so a wrong password does not trip the portal's lockout policy.
login_allowed() {
    local n ts
    n=$(sed -n 's/^fails=//p' "$STATE_FILE" 2>/dev/null | tail -1)
    ts=$(sed -n 's/^last=//p' "$STATE_FILE" 2>/dev/null | tail -1)
    if (( ${n:-0} >= 3 )) && (( $(date +%s) - ${ts:-0} < FAIL_PAUSE )); then
        log "login paused: $n consecutive failures, retry after cooldown"
        return 1
    fi
    return 0
}

write_state() {     # $1=fails $2=state — single writer for STATE_FILE
    printf 'fails=%s\nlast=%s\nstate=%s\n' "$1" "$(date +%s)" "$2" > "$STATE_FILE"
}

record_result() {   # $1 = 0 success, 1 failure; preserves state= key
    local n state
    state=$(sed -n 's/^state=//p' "$STATE_FILE" 2>/dev/null | tail -1)
    if [[ "$1" == "0" ]]; then
        n=0
    else
        n=$(sed -n 's/^fails=//p' "$STATE_FILE" 2>/dev/null | tail -1)
        n=$((${n:-0}+1))
    fi
    write_state "$n" "$state"
}

# ---------- Panabit gateway ----------
# Password encryption as implemented by the portal's own crypto.js
# (AES-128-ECB, ZeroPadding, hex output). Verified against the real JS
# via --selftest.
pa_aes_encode() {
    local text="$1" len pad
    if [ -z "$text" ]; then printf ''; return; fi
    len=$(printf '%s' "$text" | wc -c)
    pad=$(( (16 - len % 16) % 16 ))
    { printf '%s' "$text"; [ "$pad" -gt 0 ] && head -c "$pad" /dev/zero; } \
      | openssl enc -aes-128-ecb -K "$PA_KEY_HEX" -nosalt -nopad 2>/dev/null \
      | od -An -v -tx1 | tr -d ' \n'
}

pa_selftest() {
    local t w got bad=0
    while IFS='|' read -r t w; do
        got=$(pa_aes_encode "$t")
        if [[ "$got" == "$w" ]]; then
            log "selftest PASS: $(printf '%q' "$t")"
        else
            log "selftest FAIL: $(printf '%q' "$t") got=$got want=$w"
            bad=1
        fi
    done <<'VECTORS'
|
a|de285293616f1c8a2076a0cd414f8a96
abc123|3d11e715432e1241bc898117ac10f436
hello world!|bf669ad46911133b8716a0172705a5f7
sixteen chars ok|29b71c6f03826380e01ee1099747e63b
seventeen chars!!|8897d863a611710cdc3ba275e494c1cb4d3833d6591c222802fd6706f48abc5a
VECTORS
    return $bad
}

detect_portal() {
    # Sets GATEWAY (panabit|drcom) and PORTAL_BASE. When unauthenticated the
    # intercept redirect names the portal; when already online there is no
    # redirect, so fall back to PA_PORTAL (force-login case).
    local hdr loc
    hdr=$(curl -s --noproxy '*' --max-time 8 -D- -o /dev/null 'http://www.baidu.com' 2>/dev/null)
    loc=$(grep -i '^Location:' <<<"$hdr" | head -1 | tr -d '\r' | awk '{print $2}')
    if grep -qi 'panabit' <<<"$loc"; then
        GATEWAY="panabit"
        [[ "$loc" =~ ^(https?://[^/?]+) ]] && PORTAL_BASE="${BASH_REMATCH[1]}" || PORTAL_BASE="$PA_PORTAL"
    elif grep -qiE 'eportal|netauth' <<<"$hdr"; then
        GATEWAY="drcom"
        PORTAL_BASE="$EPORTAL"
    else
        GATEWAY="panabit"          # no intercept = online; default portal
        PORTAL_BASE="$PA_PORTAL"
    fi
}

panabit_login() {
    resolve_creds || return 1
    local ip mac auth_type resp rcode rmsg
    ip=$(my_ip)
    mac=$(my_mac)
    if [[ -z "$ip" || -z "$mac" || -z "$PORTAL_BASE" ]]; then
        log "ERROR: cannot determine ip/mac/portal (ip=$ip mac=$mac portal=$PORTAL_BASE)"
        return 1
    fi

    # auth_type comes from the gateway's own policy (falls back to "panabit")
    auth_type=$(curl -s --noproxy '*' --max-time 8 \
        "${PORTAL_BASE}/api?route=portal&action=load_portal_conf&ip=${ip}&mac=${mac}&device=pc" 2>/dev/null \
        | python3 -c 'import json,sys
try: print(json.load(sys.stdin).get("data",{}).get("policy",{}).get("auth1") or "panabit")
except Exception: print("panabit")' 2>/dev/null)
    auth_type="${auth_type:-panabit}"

    resp=$(curl -s --noproxy '*' --max-time 10 \
        -H 'Content-Type: application/x-www-form-urlencoded; charset=GB2312' \
        --data "u8name=${CRED_USER}" \
        "${PORTAL_BASE}/api?route=webauth&action=user_login&auth_type=${auth_type}&ip=${ip}&mac=${mac}&code=$(pa_aes_encode '')&username=${CRED_USER}&password=$(pa_aes_encode "$CRED_PASS")&remember_me=0" 2>/dev/null)

    # portal responds in GB2312; convert once to UTF-8, then parse
    local resp_u
    resp_u=$(printf '%s' "$resp" | iconv -f GB2312 -t UTF-8 -c 2>/dev/null)
    rcode=$(printf '%s' "$resp_u" | python3 -c 'import json,sys
try: print(json.load(sys.stdin).get("code",""))
except Exception: print("")' 2>/dev/null)
    rmsg=$(printf '%s' "$resp_u" | python3 -c 'import json,sys
try: print(json.load(sys.stdin).get("msg",""))
except Exception: print("")' 2>/dev/null)

    if [[ "$rcode" == "0" ]]; then
        log "Login succeed (account=${CRED_USER}, ip=${ip}, gateway=panabit)"
        return 0
    elif grep -qE '已经登录|已在线|在线' <<<"$rmsg"; then
        log "Already online (account=${CRED_USER}, ip=${ip})"
        return 0
    else
        log "Login FAILED (gateway=panabit, code=${rcode:-?}, msg=${rmsg:-<unparseable>}, resp=${resp:0:200})"
        return 1
    fi
}

# ---------- Dr.COM gateway (unchanged protocol) ----------
drcom_login() {
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
        log "Login succeed (account=${CRED_USER}, ip=${ip}, gateway=drcom)"
        return 0
    elif grep -q 'Dr.COMWebLoginID_2.htm' <<<"$resp"; then
        log "Login rejected by portal (bad credentials, or already online — check with --status)"
        return 1
    else
        log "Login FAILED: unexpected response: ${resp:0:200}"
        return 1
    fi
}

drcom_logout() {
    curl -sk --noproxy '*' --max-time 10 \
        "${EPORTAL}/?c=ACSetting&a=Logout&url=drappall" -o /dev/null
    log "Logout sent"
}

do_login() {
    detect_portal
    case "$GATEWAY" in
        panabit) panabit_login ;;
        drcom)   drcom_login ;;
    esac
}

do_logout() {
    detect_portal
    case "$GATEWAY" in
        drcom)   drcom_logout ;;
        panabit) log "Panabit portal has no logout API; nothing sent." ;;
    esac
}

# ---- argument parsing (options and command may appear in any order) ----
COMMAND="check"
ARG_USER=""; ARG_PASS=""
while (($#)); do
    case "$1" in
        --status|status)     COMMAND="status" ;;
        --login|login)       COMMAND="login" ;;
        --logout|logout)     COMMAND="logout" ;;
        --selftest|selftest) COMMAND="selftest" ;;
        check)               COMMAND="check" ;;
        -u|--user) ARG_USER="${2:?--user needs a value}"; shift ;;
        -p|--pass) ARG_PASS="${2:?--pass needs a value}"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log "Unknown argument: $1"; usage; exit 2 ;;
    esac
    shift
done

case "$COMMAND" in
    selftest)
        pa_selftest
        ;;
    status)
        if online; then log "ONLINE"; else log "OFFLINE (gateway intercepting)"; fi
        ;;
    logout)
        do_logout
        ;;
    login)
        do_login; record_result $?
        ;;
    check)
        st="offline"; online && st="online"
        prev=$(sed -n 's/^state=//p' "$STATE_FILE" 2>/dev/null | tail -1)
        if [[ "$st" != "$prev" ]]; then
            log "state: ${prev:-unknown} -> $st"
            f=$(sed -n 's/^fails=//p' "$STATE_FILE" 2>/dev/null | tail -1)
            write_state "${f:-0}" "$st"
        fi
        if [[ "$st" == "offline" ]]; then
            if ! login_allowed; then
                :   # cooldown: already logged by login_allowed
            elif do_login; then
                record_result 0
                online && log "verified ONLINE after login" || log "login ok but probe still failing"
            else
                record_result 1
            fi
        fi
        ;;
esac
