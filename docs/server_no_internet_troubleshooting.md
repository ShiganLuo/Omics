# 服务器无法访问外网：原理、排查与解决

> 2026-10-30 整理。涉及两段故障：2026-10-08 首次发现并解决"外网被劫持"；
> 2026-10-30 服务器宕机重启后旧方案失效，重新定位并修复。
> 保活脚本 `campus_auth.sh` 全文见附录。

---

## 1. 现象

服务器（IP 10.9.22.26）上访问外网出现两类异常：

- **HTTP/HTTPS 被劫持**：访问 `www.baidu.com` 等任意站点都被重定向到一个登录页；
  HTTPS 不加 `-k` 直接报证书错误——证书不是 GlobalSign 签发的真证书，而是一张伪造的
  "Web secure CA" 证书（中间人替换）。
- **DNS 伪造**：`getent hosts` 解析出的结果指向网关侧的劫持服务器，而不是真实站点。

而 SSH、内网访问一切正常（SSH 故障是另一次独立事件：10-28 深夜服务器宕机，
详见 §6 时间线，与外网问题无关）。

## 2. 原理：校园网 portal（Dr.COM）准入认证

南开校园网使用**城市热点 Dr.COM** 准入网关。其工作方式是典型的"先认证后放行"：

```
服务器 ──► 准入网关(wlanacname=blt2_) ──► 互联网
              │
              ├─ 该 IP 未认证 ──► 拦截全部 HTTP/HTTPS：
              │     · HTTP   302 重定向到 portal 登录页
              │     · HTTPS  无法解密, 用伪造证书做 MITM,
              │              假装是目标站点(故报证书错误)
              │     · DNS    伪造应答, 把域名指到劫持服务器
              │
              └─ 该 IP 已认证 ──► 透明放行, 证书/DNS 恢复正常
```

关键点：

1. **拦截不看内容只看 IP**。网关只为 IP 维护一张"已认证表"，不关心是谁在用。
2. **劫持是网络层的中间人行为**。未认证状态下 TLS 证书必然是假的——这不是站点
   出了问题，恰恰是"外网被拦截"最可靠的判别信号：假证书 + DNS 伪造同时出现，
   基本可断定是准入网关，而非路由或防火墙故障。
3. **认证 = 向 portal 发一次登录请求**，portal 校验后把该 IP 加入放行表。认证后
   HTTPS 证书立即恢复为 GlobalSign 等真实签发（10-08 实测：登录后 baidu 证书
   issuer 变回 GlobalSign RSA OV SSL CA 2018）。

### portal 登录接口（协议细节）

portal 页面本身由 `DrcomServer1.0` 提供，登录协议的参数写在门户页内嵌 JS
（`a41.js`）里，逆向后得到纯命令行可调用的接口：

```
POST http://netauth.nankai.edu.cn:801/eportal/?c=ACSetting&a=Login&url=drappall

DDDDD=,0,<账号>          # 账号节点, ",0," 是账号类型前缀(校园用户)
upass=<密码>             # 密码节点
0MKKey=123456789         # 固定值(登录按钮的伪值)
wlanacname=blt2_         # 接入控制器名
wlanuserip=<本机IP>      # 网关看到的源 IP
v6ip=                    # 留空
```

响应以 HTML 片段形式返回，标志位在注释里：

| 标志 | 含义 |
|---|---|
| `Dr.COMWebLoginID_3.htm` | 登录成功 |
| `Dr.COMWebLoginID_2.htm` | 登录失败，`msga='...'` 给出具体原因 |
| `Dr.COMWebLoginID_0.htm` | portal 登录页本身（未认证时被重定向到它） |

失败原因样例：`msga='账号不存在'`、`msga='统一身份认证验证失败，如忘记密码可通过
访问 https://iam.nankai.edu.cn ...'`、`msga='请尝试在自动弹出或访问除netauth以外
域名重定向的认证页面登录'`（这条是"已在线"时重复登录的例行拒绝）。

注销接口对称：`?c=ACSetting&a=Logout&url=drappall`（异步生效，发出后不会立即掉线）。

### portal 校验行为的两个阶段（重要）

- **2026-10-08（宽松期）**：IP 未认证时 portal **不校验账号密码**，任意假账号
  （如 FAKEUSER_TEST）都返回 `Login succeed` 并放行网络。当时就是靠假账号登上的。
- **2026-10-30（校验期）**：宕机重启后 IP 重新变为未认证，此时 portal **已改为
  真正校验账号**：假账号返回"拒绝"，连不存在的账号名直接报"账号不存在"。
  说明 portal 侧做过策略/账号库升级，免校验窗口不再存在。

**教训**：假账号能登上去只是利用了 portal 的临时宽松行为，不可依赖；
`~/.campus_auth` 里必须放能被 portal 认可的真实账号（学号/工号）。

## 3. 排查过程（方法论）

1. **判别拦截层**：`curl -s https://www.baidu.com`（不加 `-k`）看证书 issuer。
   伪造证书 + DNS 伪造 → 准入网关劫持，不是路由/防火墙/hosts 问题。
2. **找到 portal**：任意 HTTP 请求会被 302 到登录页；页面 `Server: DrcomServer1.0`、
   嵌入 JS 里有 `authloginpath='/eportal/?c=ACSetting&a=Login'` 等配置，
   据此还原登录接口（见 §2）。
3. **确认成功判据**：不要只看接口返回——登录后实测三点：
   直连 baidu HTTP 200、证书恢复真实签发、DNS 不再伪造。
4. **失败响应要看 `msga` 原文**，区分"账号不存在 / 密码错 / 已在线"三种情况，
   否则会像本次一样把"账号不存在"误判成"密码不对"。

## 4. 解决方案

### 4.1 登录脚本 `campus_auth.sh`

单文件 bash 脚本（全文见附录），功能：探测外网连通性 → 掉线时自动向 portal
重登 → 凭据持久化。

```bash
campus_auth.sh                  # 默认模式：在线不动，掉线自动登录（cron 用）
campus_auth.sh --status         # 只查状态，绝不登录
campus_auth.sh --login          # 强制登录（用 ~/.campus_auth 里的凭据）
campus_auth.sh --logout         # portal 注销
campus_auth.sh --login -u <账号> -p <密码>   # 登录并把凭据存入 ~/.campus_auth
```

实现要点：

- 在线探测：`curl https://www.baidu.com` 返回 200 即认为已认证（劫持状态下
  不会是 200）。
- 登录请求即 §2 的 POST；响应含 `ID_3` = 成功，含 `ID_2` = 失败（进一步读 `msga`）。
- 凭据优先级：`-u/-p` > `~/.campus_auth`（chmod 600）> 交互式隐藏输入。
  传入的 `-u/-p` 会自动写入 `~/.campus_auth`，供 cron 复用。
- cron 空环境下 HOME 可能未定义，脚本内 `getent passwd` 兜底。

### 4.2 cron 保活

```cron
*/10 * * * * <脚本绝对路径>/campus_auth.sh >> $HOME/.campus_auth.log 2>&1
```

日志 `~/.campus_auth.log`。正常时每 10 分钟一行 `ONLINE`；
出现 `OFFLINE, attempting portal login...` 即发生过掉线并触发自动重登。

### 4.3 本次（10-30）的最终解决

```
campus_auth.sh --login -u <学号> -p <真实密码>
→ portal 返回 ID_3, 网络放行（baidu 200）
→ 凭据写入 ~/.campus_auth，cron 保活恢复自动运行
```

候选密码逐一实测：第 1 个报"统一身份认证验证失败"（密码错）；
第 2 个成功；第 3 个未轮到（IP 已在线后重复登录会被例行拒绝，不再反映密码真伪）。

## 5. 操作手册（下次掉线照做）

1. `campus_auth.sh --status` —— 确认 OFFLINE。
2. `campus_auth.sh --login` —— 用已存凭据重登；若报"账号不存在"才需要
   `-u <真实账号> -p <密码>` 重新指定。
3. `campus_auth.sh --status` + `curl -s https://www.baidu.com` 验证放行。
4. 若反复失败，读 `msga` 原文区分：账号不存在（账号名/账号状态问题，找网络中心，
   统一身份认证入口 https://iam.nankai.edu.cn 可自助找回密码）/ 密码错 / 已在线
   （正常现象，稍后 --status 复核即可）。
5. 不要再用假账号试探——10-30 起 portal 会真正校验。

## 6. 时间线（含 10-28 宕机，供对照）

| 时间 | 事件 |
|---|---|
| 10-08 21:08 | 发现 HTTPS 劫持（假证书 + DNS 伪造），逆向 Dr.COM portal 登录接口 |
| 10-08 21:08 | 用假账号登录成功（当时 portal 不校验），网络放行，证书恢复 |
| 10-08 21:12 | 写出 campus_auth.sh，装 cron 每 10 分钟保活 |
| 10-08 ~ 10-28 | 持续 ONLINE（`~/.campus_auth.log` 每 10 分钟一条），仅 10-14 一次瞬时掉线自愈 |
| 10-28 22:40~22:50 | **服务器宕机**（last 记录为 crash，非正常关机；宕机原因需 root 看内核日志），10-29 全天不可达（SSH 连不上的原因） |
| 10-30 07:51 | 服务器重新开机 |
| 10-30 08:10 起 | cron 自动重登持续失败：portal 转为校验账号，假账号被拒、luosg 报"账号不存在" |
| 10-30 11:26 | 用真实学号账号登录成功，网络恢复，凭据持久化 |

## 7. 经验总结

1. **假证书 + DNS 伪造 = 准入网关劫持**，不是站点故障；先看证书 issuer 再查路由。
2. portal 的免校验窗口是临时行为，**凭据必须是 portal 认可的真实账号**（学号/工号）。
3. 登录失败必须读 `msga` 原文再下结论——"账号不存在"曾被脚本笼统归为
   "bad credentials, or already online"，误导排查方向。
4. 保活脚本 + cron 是长效方案；日志（`~/.campus_auth.log`）是掉线历史的第一证据。
5. `--logout` 异步生效，注销后立刻 `--status` 仍显示 ONLINE 是正常现象。

---

## 附录：campus_auth.sh 全文

```bash
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
#   */10 * * * * <脚本路径>/campus_auth.sh >> $HOME/.campus_auth.log 2>&1

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
```
