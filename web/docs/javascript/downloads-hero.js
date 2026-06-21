(function() {
'use strict';

function initHeroDownload()
{
  const el = document.getElementById('hero-download');
  if (!el || !window.PythonSCADDownloads) {
    return;
  }

  const {
    fetchLatestRelease,
    groupAssetsByPlatform,
    detectUserPlatform,
    pickAssetForPlatform,
    PLATFORM_LABELS
  } = window.PythonSCADDownloads;

  el.textContent = 'Loading latest release…';

  async function renderHeroDownload()
  {
    try {
      const data = await fetchLatestRelease();
      const {byPlatform} = groupAssetsByPlatform(data.assets || []);
      const userPlatform = detectUserPlatform();

      let html =
        `<p class="hero-download-version">Latest version: <strong>${data.tag_name}</strong></p>`;
      html += `<p class="hero-download-actions">`;

      if (userPlatform === 'linux-debian') {
        html += `<a class="md-button md-button--primary hero-download-button" `;
        html += `href="installation/#linux-debian-ubuntu-apt">Install via APT</a>`;
        html += ` <a class="md-button hero-download-all" href="downloads/">All downloads</a>`;
        html += `</p>`;
        html += `<p class="hero-download-file">Recommended: `;
        html += `<a href="https://repos.pythonscad.org/apt/">repos.pythonscad.org/apt</a></p>`;
      } else if (userPlatform === 'linux-fedora') {
        html += `<a class="md-button md-button--primary hero-download-button" `;
        html += `href="installation/#linux-fedora-centos-rhel-yum-dnf">Install via YUM/DNF</a>`;
        html += ` <a class="md-button hero-download-all" href="downloads/">All downloads</a>`;
        html += `</p>`;
        html += `<p class="hero-download-file">Recommended: `;
        html += `<a href="https://repos.pythonscad.org/yum/">repos.pythonscad.org/yum</a></p>`;
      } else {
        const picked = pickAssetForPlatform(byPlatform, userPlatform);

        if (picked) {
          const label = PLATFORM_LABELS[picked.platform] || picked.platform;
          html += `<a class="md-button md-button--primary hero-download-button" `;
          html += `href="${picked.asset.browser_download_url}">`;
          html += `Download for ${label}</a>`;
          html += ` <a class="md-button hero-download-all" href="downloads/">All downloads</a>`;
          html += `</p>`;
          html += `<p class="hero-download-file">${picked.asset.name}</p>`;
        } else {
          html += `<a class="md-button md-button--primary hero-download-button" href="downloads/">`;
          html += `Download PythonSCAD</a>`;
          html += `</p>`;
        }
      }

      el.innerHTML = html;
    } catch (e) {
      el.innerHTML = `<p class="hero-download-error">Could not load release info. ` +
        `<a href="downloads/">Browse downloads</a>.</p>`;
    }
  }

  renderHeroDownload();
}

if (window.PythonSCADDownloads) {
  window.PythonSCADDownloads.onMkDocsPageLoad(initHeroDownload);
} else {
  document.addEventListener('DOMContentLoaded', initHeroDownload);
}
})();
