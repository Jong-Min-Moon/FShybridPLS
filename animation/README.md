# Hybrid PLS animation

Educational walkthrough of Hybrid PLS (motivation → associations → hybrid \(W\) → PLS weights → one score \(\rho\)).

## View

- Open [`index.html`](index.html) in a browser (loads sibling `animation_data.js`).
- README on GitHub shows the looping demo GIF [`hybridpls-demo.gif`](hybridpls-demo.gif) (GitHub cannot run interactive HTML inside README.md).

## Record / refresh the GIF

```bash
cd animation
npm install
npx playwright install chromium
npm run record-gif
```

This writes `hybridpls-demo.gif` (~3 MB). Use `index.html?record=1` to hide page chrome while capturing.

## GitHub Pages

**One-time setup (required before the workflow can succeed):**

1. Open **Settings → Pages** for this repo:  
   https://github.com/Jong-Min-Moon/FShybridPLS/settings/pages
2. Under **Build and deployment → Source**, choose **GitHub Actions** (not “Deploy from a branch”).
3. Re-run the failed **pages-animation** workflow (Actions → pages-animation → Re-run jobs), or push any change under `animation/`.

After that, the interactive demo is at `https://jong-min-moon.github.io/FShybridPLS/`.

Pushing changes under `animation/` runs `.github/workflows/pages-animation.yml`. The first failure with `Get Pages site failed … Not Found` almost always means step 2 above was not done yet — `GITHUB_TOKEN` cannot enable Pages by itself.