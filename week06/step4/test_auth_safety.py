"""Offline security checks with synthetic credentials; no network or real login."""
import contextlib
import io
import json
import unittest
from unittest.mock import patch
import urllib.error
import prepare_cdse as target


class AuthSafetyTests(unittest.TestCase):
    def execute_fake_auth(self, transport):
        output = io.StringIO()
        with patch.object(target.sys.stdin, "isatty", return_value=True), \
             patch.object(target.sys.stdout, "isatty", return_value=True), \
             patch.object(target.getpass, "getpass", side_effect=["SYNTHETIC_EMAIL", "SYNTHETIC_PASSWORD", "123456"]), \
             patch.object(target, "opener", return_value=transport), contextlib.redirect_stdout(output):
            # redirect_stdout replaces the patched stream, so explicitly mark it interactive.
            output.isatty = lambda: True
            result = target.check_authentication()
        return result, output.getvalue()

    def test_success_uses_only_token_endpoint_and_never_prints_secrets(self):
        class FakeTransport:
            def open(self, request, timeout):
                self.request = request
                assert request.full_url == target.TOKEN_URL
                assert request.get_method() == "POST"
                assert b"totp=123456" in request.data
                return contextlib.closing(io.BytesIO(json.dumps({
                    "access_token": "SYNTHETIC_ACCESS_TOKEN", "refresh_token": "SYNTHETIC_REFRESH_TOKEN", "expires_in": 600}).encode()))
        transport = FakeTransport()
        status, output = self.execute_fake_auth(transport)
        self.assertEqual(status, 0)
        for secret in ["SYNTHETIC_EMAIL", "SYNTHETIC_PASSWORD", "SYNTHETIC_ACCESS_TOKEN", "SYNTHETIC_REFRESH_TOKEN", "123456"]:
            self.assertNotIn(secret, output)
        self.assertIsNone(transport.request.data)

    def test_http_error_does_not_expose_body_or_headers(self):
        class FakeTransport:
            def open(self, request, timeout):
                raise urllib.error.HTTPError(target.TOKEN_URL, 401,
                    "SYNTHETIC_PASSWORD", {"Secret": "SYNTHETIC_ACCESS_TOKEN"}, io.BytesIO(b"SYNTHETIC_PASSWORD"))
        with self.assertRaises(target.PreparationError) as failure:
            self.execute_fake_auth(FakeTransport())
        self.assertNotIn("SYNTHETIC", str(failure.exception))
        self.assertIn("HTTP 401", str(failure.exception))

    def test_noninteractive_auth_is_refused_before_input_or_network(self):
        with patch.object(target.sys.stdin, "isatty", return_value=False), \
             patch.object(target.getpass, "getpass") as prompt, patch.object(target, "opener") as network:
            with self.assertRaises(target.PreparationError):
                target.check_authentication()
            prompt.assert_not_called()
            network.assert_not_called()

    def test_auth_redirects_are_not_followed(self):
        self.assertIsNone(target.NoRedirect().redirect_request(None, None, 302, "", {}, "https://example.invalid"))

    def test_public_preflight_cannot_call_download_or_token_endpoint(self):
        with patch.object(target, "opener") as network:
            for url in [target.TOKEN_URL, "https://download.dataspace.copernicus.eu/odata/v1/Products(x)/$value"]:
                with self.assertRaises(target.PreparationError):
                    target.public_json(url)
            network.assert_not_called()


if __name__ == "__main__":
    unittest.main(verbosity=2)
