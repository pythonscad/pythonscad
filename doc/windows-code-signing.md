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
the **public certificate** (PEM format) for the CI workflow.

### 2.1 Install proCertum CardManager

Download and install **proCertum CardManager** from the Certum website. This
provides the PKCS#11 module that jsign uses to communicate with SimplySign.

### 2.2 Export the certificate as PEM

On Windows, open a PowerShell prompt and run:

```powershell
# List certificates on your SimplySign account
certutil -scinfo

# Export the certificate to PEM (replace CERT_THUMBPRINT with your cert's thumbprint)
$cert = Get-ChildItem Cert:\CurrentUser\My | Where-Object { $_.Thumbprint -eq 'CERT_THUMBPRINT' }
$pemContent = "-----BEGIN CERTIFICATE-----`n" +
  [Convert]::ToBase64String($cert.RawData, 'InsertLineBreaks') +
  "`n-----END CERTIFICATE-----"
$pemContent | Out-File -FilePath certum.pem -Encoding ascii
```

Alternatively, use OpenSSL once the certificate appears in your system store:

```bash
openssl x509 -in certum.cer -inform DER -out certum.pem
```

### 2.3 Base64-encode the PEM for GitHub

```bash
base64 -w 0 certum.pem
```

Copy the entire output — you will paste it into a GitHub secret.

### 2.4 Find your certificate's CN

The MSIX `publisher-cn` field must exactly match the Subject CN of your
certificate. Find it with:

```powershell
$cert = Get-ChildItem Cert:\CurrentUser\My | Where-Object { $_.Thumbprint -eq 'CERT_THUMBPRINT' }
$cert.Subject
```

It will look something like `CN=Your Name, O=Open Source Developer, ...`.
You need the full Subject string in exactly this format.

---

## Part 3: GitHub Secrets

In your repository: **Settings → Secrets and variables → Actions**, add:

| Secret name | Description |
| --- | --- |
| `CERTUM_CERTIFICATE_BASE64` | Full base64 output of your `certum.pem` file (from step 2.3) |
| `CERTUM_CARD_PIN` | Your SimplySign card PIN |
| `CERTUM_CERTIFICATE_CN` | Full Subject string of your certificate (from step 2.4), used as MSIX `publisher-cn` |

The workflow skips signing silently when `CERTUM_CERTIFICATE_BASE64` or
`CERTUM_CARD_PIN` are absent, so unsigned test builds continue to work without
any changes.

---

## Part 4: How the workflow uses this

The `Build Windows Native Package` workflow
(`.github/workflows/build-windows-native.yml`) calls the reusable composite
action `.github/actions/sign-windows/action.yml` twice:

1. **After NSIS packaging** — signs `artifacts/*.exe` (the installer)
2. **After MSIX packaging** — signs the `.msix` file

When `CERTUM_CERTIFICATE_CN` is also set, the MSIX `publisher-cn` is
automatically switched from the unsigned-namespace OID placeholder to your real
certificate subject, which is required for a validly signed MSIX.

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
3. Install proCertum CardManager; export your public certificate as PEM.
4. Base64-encode the PEM and add it as `CERTUM_CERTIFICATE_BASE64` secret.
5. Add your SimplySign PIN as `CERTUM_CARD_PIN` secret.
6. Add your certificate Subject string as `CERTUM_CERTIFICATE_CN` secret.
7. Trigger the Windows package build; approve the SimplySign notification(s).
8. Verify signatures on the downloaded artifacts.

Keep your SimplySign PIN secret. Never commit certificate private keys — with
SimplySign there are none to commit.

<!-- markdownlint-enable MD013 -->
