# Windows Code Signing with Certum SimplySign

<!-- markdownlint-disable MD013 -->

This guide explains how to sign PythonSCAD Windows builds (NSIS installer `.exe`
and MSIX package) using a **Certum Open Source Code Signing** certificate on
the SimplySign cloud HSM. The private key never leaves Certum's infrastructure;
signing is authorized via a push notification on your Android or iOS device.

## Overview

- Signing uses [jsign](https://ebourg.github.io/jsign/), a Java-based
  Authenticode signing tool with native Certum SimplySign support.
- The GitHub Actions workflow signs artifacts automatically during the
  `Build Windows Native Package` job.
- When a signing operation is triggered, the **SimplySign mobile app** sends a
  push notification to your phone. Approving it (with your PIN or biometric)
  unblocks the workflow and completes the signature.
- If the signing secrets are not configured the workflow still runs and produces
  unsigned artifacts — suitable for snapshot/test builds.

## Part 1: Purchase the Certificate

Buy the **Open Source Code Signing on SimplySign** product from
[certum.store](https://certum.store/open-source-code-signing-on-simplysign.html).

- ~$58/year (pricing may change)
- 5 000 signatures per month
- Requires proof of open-source project status

During the purchase/onboarding process Certum will:

1. Verify your identity and project
2. Issue a certificate tied to your SimplySign account
3. Guide you through installing the **SimplySign** app on your Android or iOS
   device

## Part 2: Export your public certificate

The private key stays on Certum's HSM and cannot be exported. You only need
the **public certificate in DER format** (base64-encoded) for the CI workflow.

### 2.1 Install proCertum CardManager

Download and install **proCertum CardManager** from the Certum website. This
provides the PKCS#11 module that jsign uses to communicate with SimplySign.

### 2.2 Export the certificate as DER

The CI action expects **base64-encoded DER** (the raw binary certificate format),
not a PEM text file. Export the DER bytes from your system certificate store:

On Windows (PowerShell):

```powershell
# List certificates — find your Certum cert by thumbprint or issuer
certutil -scinfo

# Export raw DER bytes (replace CERT_THUMBPRINT with your cert's thumbprint)
$cert = Get-ChildItem Cert:\CurrentUser\My | Where-Object { $_.Thumbprint -eq 'CERT_THUMBPRINT' }
[IO.File]::WriteAllBytes("certum.cer", $cert.RawData)
```

On Linux/macOS (if the cert was already downloaded as `.pem`):

```bash
openssl x509 -in certum.pem -outform DER -out certum.cer
```

### 2.3 Base64-encode the DER for GitHub

```bash
# Linux (GNU coreutils)
base64 -w0 certum.cer

# macOS (BSD base64 — does not support -w)
base64 certum.cer | tr -d '\n'

# PowerShell (Windows)
[Convert]::ToBase64String([IO.File]::ReadAllBytes("certum.cer"))
```

Copy the entire single-line output — you will paste it into a GitHub secret.

**Important:** the secret must be base64 of the raw DER bytes, not base64 of a
PEM text file. The two look similar but are not interchangeable.

### 2.4 Find your certificate's full Subject DN

The MSIX `publisher-cn` field must exactly match the full **Subject Distinguished
Name (DN)** of your certificate — not just the CN component. Find it with:

```powershell
$cert = Get-ChildItem Cert:\CurrentUser\My | Where-Object { $_.Thumbprint -eq 'CERT_THUMBPRINT' }
$cert.Subject
```

It will look something like `CN=Your Name, O=Open Source Developer, C=PL`.
Copy the **entire string** exactly as printed — all components, in the same order.

---

## Part 3: GitHub Secrets

In your repository: **Settings → Secrets and variables → Actions**, add:

| Secret name | Description |
| --- | --- |
| `CERTUM_CERTIFICATE_BASE64` | Full base64 output of your `certum.cer` DER file (from step 2.3) |
| `CERTUM_CARD_PIN` | Your SimplySign card PIN |
| `CERTUM_CERTIFICATE_CN` | Full Subject DN of your certificate (from step 2.4), used as MSIX `publisher-cn` |

The workflow logs a skip message and produces unsigned artifacts when any of
the required secrets are absent, so test builds continue to work without
any changes.

---

## Part 4: How the workflow uses this

The `Build Windows Native Package` workflow
(`.github/workflows/build-windows-native.yml`) calls the reusable composite
action `.github/actions/sign-windows/action.yml` twice:

1. **After NSIS packaging** — signs the resolved installer path
2. **After MSIX packaging** — signs the `.msix` file

When all three secrets (`CERTUM_CERTIFICATE_BASE64`, `CERTUM_CARD_PIN`, and
`CERTUM_CERTIFICATE_CN`) are present, the MSIX `Identity Publisher` field is
automatically set to your real certificate Subject DN. Otherwise the
unsigned-namespace OID placeholder is used, which is required for a validly
signed MSIX — the DN must match the signing certificate exactly.

### Runner prerequisite: proCertum CardManager

`jsign --storetype CRYPTOCERTUM` communicates with the SimplySign cloud HSM
via the **proCertum CardManager PKCS#11 driver** (`sc30pkcs11.dll`).
GitHub-hosted `windows-latest` runners do not include this software.

The sign-windows action checks for the DLL at startup and fails fast with a
clear error if it is missing. Before enabling signing in production you must
add a silent-install step to `.github/actions/sign-windows/action.yml`:

```powershell
Invoke-WebRequest <pinned-url> -OutFile cardmanager.exe
if ((Get-FileHash cardmanager.exe -Algorithm SHA256).Hash.ToLower() -ne '<sha256>') { exit 1 }
Start-Process cardmanager.exe -ArgumentList '/S' -Wait
```

Download the installer from the Certum website, pin both the URL and SHA-256,
and replace the placeholder comment in the action before your first signed release.

### Mobile app interaction

Each file signing operation contacts the SimplySign cloud HSM. If your account
is configured for push-notification authorization:

- The **SimplySign app** shows a notification on your phone
- You approve with your PIN or biometric
- The workflow step unblocks and continues

For a release with two signed files (`.exe` and `.msix`) you will receive
**two notifications**. Keep the GitHub Actions page open in a browser tab so
you can see progress.

### Timeout

GitHub Actions steps time out after **6 hours** by default — you have plenty
of time to respond to the phone notification at your convenience during a
release.

---

## Part 5: Verifying signatures

After downloading a release artifact, verify the signature on Windows:

```powershell
# Check the NSIS installer
Get-AuthenticodeSignature .\PythonSCAD-*-Installer.exe | Format-List

# Check the MSIX
Get-AuthenticodeSignature .\PythonSCAD-*.msix | Format-List
```

Both should show `Status: Valid` and a `SignerCertificate` issued by Certum.

On Linux with osslsigncode:

```bash
osslsigncode verify -in PythonSCAD-*-Installer.exe
```

---

## Summary checklist

1. Purchase Certum Open Source Code Signing on SimplySign.
2. Install SimplySign app on Android or iOS and link it to your account.
3. Install proCertum CardManager; export your public certificate as DER (`certum.cer`).
4. Base64-encode the DER file and add it as `CERTUM_CERTIFICATE_BASE64` secret.
5. Add your SimplySign PIN as `CERTUM_CARD_PIN` secret.
6. Add your certificate Subject string as `CERTUM_CERTIFICATE_CN` secret.
7. Trigger the Windows package build; approve the SimplySign notification(s).
8. Verify signatures on the downloaded artifacts.

Keep your SimplySign PIN secret. Never commit certificate private keys — with
SimplySign there are none to commit.

<!-- markdownlint-enable MD013 -->
