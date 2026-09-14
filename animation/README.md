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

Pushing changes under `animation/` runs `.github/workflows/pages-animation.yml`. Enable **Settings → Pages → Source: GitHub Actions** once; then the interactive demo is at `https://jong-min-moon.github.io/FShybridPLS/`.
