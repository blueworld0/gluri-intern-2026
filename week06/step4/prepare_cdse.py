"""CDSE preparation only: local/public preflight or local hidden-input auth check.

There is deliberately NO image download implementation. No secrets in CLI args,
environment variables, files, response logs, exception messages, or stdout.
Uses only the Python standard library; run directly in a user's own terminal.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import getpass
import hashlib
import json
import platform
import shutil
import ssl
import sys
import urllib.error
import urllib.parse
import urllib.request
import uuid
import warnings

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRODUCT_ID = "6c4173d8-59e8-41f8-b339-c49620900f8f"
PRODUCT_NAME = "S2C_MSIL2A_20250428T020711_N0511_R103_T52SDD_20250428T054014.SAFE"
TOKEN_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
CATALOGUE_URL = f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products({PRODUCT_ID})?%24expand=Attributes"
DISCOVERY_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/.well-known/openid-configuration"


class PreparationError(Exception):
    """Only fixed safe messages may be placed here."""


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None  # In particular, never forward a password POST elsewhere.


def opener():
    return urllib.request.build_opener(
        urllib.request.HTTPSHandler(context=ssl.create_default_context()), NoRedirect())


def public_json(url):
    # Fixed allowlist prevents accidental calls to download/$value or auth POST.
    if url not in {CATALOGUE_URL, DISCOVERY_URL}:
        raise PreparationError("Unexpected public metadata URL; stopped.")
    try:
        request = urllib.request.Request(url, headers={"Accept": "application/json"})
        with opener().open(request, timeout=30) as response:
            payload = response.read(2_000_001)
            if len(payload) > 2_000_000:
                raise PreparationError("Public metadata exceeded size limit.")
            return json.loads(payload)
    except urllib.error.HTTPError as exc:
        # Do not include response bodies, headers or arbitrary URLs in errors.
        raise PreparationError(f"Public metadata HTTP {exc.code}; no authentication attempted.") from None
    except (urllib.error.URLError, TimeoutError, OSError, ValueError):
        raise PreparationError("Public metadata connection/format check failed; no authentication attempted.") from None


def local_preflight():
    config = json.loads((HERE / "product.json").read_text(encoding="utf-8"))
    if config["product_uuid"] != PRODUCT_ID or config["product_name"] != PRODUCT_NAME:
        raise PreparationError("Product configuration mismatch.")
    download_dir = HERE / "downloads"
    if download_dir.resolve() != HERE.resolve() / "downloads":
        raise PreparationError("Download directory redirects outside Step4; stopped.")
    download_dir.mkdir(exist_ok=True)
    disk = shutil.disk_usage(download_dir)
    zip_path = download_dir / (PRODUCT_NAME + ".zip")
    safe_path = download_dir / PRODUCT_NAME
    part_path = download_dir / (PRODUCT_NAME + ".zip.part")
    # Exclusive temporary write probe never overwrites an existing file.
    probe = download_dir / (".write-check-" + uuid.uuid4().hex)
    try:
        with probe.open("xb") as f:
            f.write(b"step4-write-check")
        probe.unlink()
        writable = True
    except OSError:
        writable = False
    old_manifest = ROOT / "week06/step3/download_manifest_primary.csv"
    with old_manifest.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    return {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "python_version": platform.python_version(),
        "python_executable": sys.executable, "openssl": ssl.OPENSSL_VERSION,
        "dependencies": "Python standard library only; no installation required",
        "download_directory": str(download_dir.resolve()),
        "future_archive_path": str(zip_path), "future_safe_path": str(safe_path),
        "directory_writable": writable, "free_space_bytes": disk.free,
        "minimum_free_space_bytes": config["minimum_free_space_bytes"],
        "has_planning_space": disk.free >= config["minimum_free_space_bytes"],
        "archive_exists": zip_path.exists(), "safe_directory_exists": safe_path.exists(),
        "partial_archive_exists": part_path.exists(),
        "step3_planned_subset_files": len(rows),
        "step3_planned_subset_bytes": sum(int(r["bytes"]) for r in rows),
        "step3_manifest_sha256": hashlib.sha256(old_manifest.read_bytes()).hexdigest(),
        "credential_input": "hidden local interactive input only",
        "authentication_tested": False, "authenticated_download_access_tested": False,
        "download_started": False, "download_implementation_present": False,
    }


def run_preflight(public=False):
    result = local_preflight()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    run_dir = HERE / "checks" / stamp
    run_dir.mkdir(parents=True, exist_ok=False)
    if public:
        for key, url in [("product", CATALOGUE_URL), ("oidc_discovery", DISCOVERY_URL)]:
            try:
                data = public_json(url)
                if key == "product":
                    if data.get("Id") != PRODUCT_ID or data.get("Name") != PRODUCT_NAME:
                        raise PreparationError("Public product identity mismatch.")
                    with (run_dir / "product_metadata.json").open("x", encoding="utf-8") as f:
                        json.dump(data, f, ensure_ascii=False, indent=2)
                    result.update(product_name=data["Name"], product_uuid=data["Id"],
                        product_online=data.get("Online"), catalogue_content_length_bytes=data.get("ContentLength"),
                        product_checksums=data.get("Checksum"), product_metadata_http_status=200)
                else:
                    if data.get("token_endpoint") != TOKEN_URL:
                        raise PreparationError("OIDC token endpoint mismatch.")
                    with (run_dir / "oidc_public_metadata.json").open("x", encoding="utf-8") as f:
                        json.dump(data, f, indent=2)
                    result.update(oidc_endpoint_verified=True,
                        oidc_password_grant_advertised="password" in data.get("grant_types_supported", []))
            except PreparationError as exc:
                result[key + "_check_error"] = str(exc)
    with (run_dir / "preflight.json").open("x", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    print("Saved public/local preflight report:", run_dir / "preflight.json")
    return 0 if (result["directory_writable"] and result["has_planning_space"]
                 and not any(k.endswith("_check_error") for k in result)) else 1


def check_authentication():
    # Never fall back to echoing a password under a pipe, notebook, or agent PTY.
    if not sys.stdin.isatty() or not sys.stdout.isatty():
        raise PreparationError("Run --check-auth directly in your own interactive PowerShell terminal; piped input is refused.")
    print("Authentication check only. No product data will be requested or downloaded.")
    request = None
    form = {}
    token_response = {}
    encoded = None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", getpass.GetPassWarning)
            form["username"] = getpass.getpass("CDSE username/email (hidden): ")
            form["password"] = getpass.getpass("CDSE password (hidden): ")
            otp = getpass.getpass("2FA current code (hidden; Enter if not enabled): ")
        if not form["username"] or not form["password"]:
            raise PreparationError("Empty local credentials; no request sent.")
        if otp:
            form["totp"] = otp
        otp = None
        form.update(client_id="cdse-public", grant_type="password")
        encoded = urllib.parse.urlencode(form).encode("utf-8")
        request = urllib.request.Request(TOKEN_URL, data=encoded, method="POST",
            headers={"Content-Type": "application/x-www-form-urlencoded", "Accept": "application/json"})
        try:
            with opener().open(request, timeout=30) as response:
                body = response.read(1_000_001)
                if len(body) > 1_000_000:
                    raise PreparationError("Authentication response too large; stopped.")
                token_response = json.loads(body)
                body = None
        except urllib.error.HTTPError as exc:
            raise PreparationError(f"Authentication failed (HTTP {exc.code}); check local credentials/2FA/account. Response body suppressed.") from None
        except (urllib.error.URLError, TimeoutError, OSError, ValueError):
            raise PreparationError("Authentication connection/format check failed. Details suppressed.") from None
        if not isinstance(token_response.get("access_token"), str) or not token_response["access_token"]:
            raise PreparationError("No access token in response; response details suppressed.")
        expires = token_response.get("expires_in")
        suffix = f" (expires_in={expires} seconds)" if isinstance(expires, int) else ""
        print("Authentication succeeded" + suffix + ". Token is not displayed or saved.")
        print("No product request was made. This process discards the token on exit; authenticate again in a future download process.")
        return 0
    finally:
        # Limit lifetimes; Python immutable bytes/strings cannot guarantee zeroization.
        form.clear()
        token_response.clear()
        encoded = None
        if request is not None:
            request.data = None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--public-preflight", action="store_true", help="GET public product and OIDC metadata only; no auth")
    mode.add_argument("--check-auth", action="store_true", help="Local hidden-input token issuance check only; no download")
    args = parser.parse_args()
    try:
        return check_authentication() if args.check_auth else run_preflight(args.public_preflight)
    except PreparationError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (EOFError, KeyboardInterrupt, getpass.GetPassWarning):
        print("Local input cancelled or hidden input unavailable; stopped.", file=sys.stderr)
        return 1
    except Exception:
        # A credential-bearing exception is never rendered or traced.
        print("Preparation stopped. Error details suppressed; no credentials or response bodies logged.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
