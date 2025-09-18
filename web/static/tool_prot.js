var currentGene = null;
let viewer = null;

// === Color Maps ===
const viridis = [
  [68, 1, 84], [71, 44, 122], [59, 81, 139], [44, 113, 142],
  [33, 144, 141], [39, 173, 129], [92, 200, 99], [170, 220, 50],
  [253, 231, 37]
];

const rdbu = [
  [5, 48, 97],
  [33, 102, 172],
  [67, 147, 195],
  [146, 197, 222],
  [247, 247, 247],
  [244, 165, 130],
  [214, 96, 77],
  [178, 24, 43],
  [103, 0, 31]
];

// === Centralized opacity lookup ===
const SURFACE_OPACITY = {
  "plddt": 0.8,
  "fpocket": 0.8,
  "sasa": 0.88,
  "pi": 0.88
};

// === Color Functions ===
function viridisColor(value) {
  if (typeof value !== "number" || isNaN(value)) value = 0.0;
  value = Math.min(1, Math.max(0, value));
  const idx = value * (viridis.length - 1);
  let i = Math.floor(idx);
  let f = idx - i;
  if (i >= viridis.length - 1) {
    i = viridis.length - 1;
    f = 0;
  }
  const next = viridis[i + 1] || viridis[i];
  const rgb = [0, 1, 2].map(j => Math.round(viridis[i][j] + f * (next[j] - viridis[i][j])));
  return (rgb[0] << 16) | (rgb[1] << 8) | rgb[2];
}

function rdbuColor(value) {
  if (typeof value !== "number" || isNaN(value)) value = 0.5;
  value = Math.min(1, Math.max(0, value));
  const idx = value * (rdbu.length - 1);
  let i = Math.floor(idx);
  let f = idx - i;
  if (i >= rdbu.length - 1) {
    i = rdbu.length - 1;
    f = 0;
  }
  const next = rdbu[i + 1] || rdbu[i];
  const rgb = [0, 1, 2].map(j => Math.round(rdbu[i][j] + f * (next[j] - rdbu[i][j])));
  return (rgb[0] << 16) | (rgb[1] << 8) | rgb[2];
}

// === Legend Drawing ===
function drawLegend(colorFunc) {
  const canvas = document.getElementById("colorLegend");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;
  if (!width || !height) return;
  ctx.clearRect(0, 0, width, height);
  for (let i = 0; i < width; i++) {
    const t = i / (width - 1);
    const rgbInt = colorFunc(t);
    ctx.fillStyle = "#" + rgbInt.toString(16).padStart(6, "0");
    ctx.fillRect(i, 0, 1, height);
  }
}

function updateLegend(gene, style) {
  if (!protData || !protData[gene]) return;
  const legend = document.getElementById("legend");
  if (!legend) return;
  const minLabel = document.getElementById("legendMin");
  const maxLabel = document.getElementById("legendMax");
  legend.style.display = "none";

  const stats = protData[gene]?.stats || null;

  switch (style) {
    case "surface-plddt":
      minLabel.innerHTML = "<strong>Low</strong> (0&#37;)";
      maxLabel.innerHTML = "<strong>High</strong> (100&#37;)";
      legend.style.display = "flex";
      drawLegend(viridisColor);
      break;

    case "surface-fpocket":
      if (stats?.fpocket) {
        minLabel.innerHTML = "<strong>Min</strong> (" + stats.fpocket.min.toFixed(3) + ")";
        maxLabel.innerHTML = "<strong>Max</strong> (" + stats.fpocket.max.toFixed(3) + ")";
        legend.style.display = "flex";
        drawLegend(viridisColor);
      }
      break;

    case "surface-sasa":
      if (stats?.sasa) {
        minLabel.innerHTML =
          "<strong>Min</strong> (" + stats.sasa.min.toFixed(2) + " &Aring;&sup2;/res)";
        maxLabel.innerHTML =
          "<strong>Max</strong> (" + stats.sasa.max.toFixed(2) + " &Aring;&sup2;/res)";
        legend.style.display = "flex";
        drawLegend(viridisColor);
      }
      break;

    case "surface-pi":
      if (stats?.pi) {
        minLabel.innerHTML = "<strong>-1.0</strong> (apolar/hydrophobic)";
        maxLabel.innerHTML = "<strong>+1.0</strong> (polar/hydrophilic)";
        legend.style.display = "flex";
        drawLegend(rdbuColor);
      }
      break;
  }
}

// === Surface Coloring ===
function colorSurface(contributions, useDiverging = false, useTransparency = false, mode = "generic") {
  if (!contributions || contributions.length === 0) return;

  const scoreMap = {};
  contributions.forEach(c => {
    if (c.residue_no == null) return;
    const numScore = Number(c.score);
    if (!Number.isFinite(numScore)) return;
    scoreMap[c.residue_no] = Math.min(Math.max(numScore, 0), 1);
  });

  const surfaceType = (mode === "sasa") ? $3Dmol.SurfaceType.SAS : $3Dmol.SurfaceType.VDW;
  const globalOpacity = SURFACE_OPACITY[mode] || 1.0;

  viewer.addSurface(surfaceType, {
    transparent: true,
    opacity: globalOpacity,
    colorfunc: function(atom) {
      if (!atom || atom.resi === undefined) return 0xAAAAAA;
      const s = scoreMap[atom.resi];
      if (typeof s !== "number" || isNaN(s)) return 0xAAAAAA;
      return useDiverging ? rdbuColor(s) : viridisColor(s);
    }
  });
}

// === Hover Labels (Theme-Aware) ===
function enableResidueHoverLabels() {
  viewer.setHoverable({}, true,
    function(atom, viewerInstance) {
      if (!atom || atom.resi === undefined) return;

      // Read CSS variables for theme colors
      const bgColor = getComputedStyle(document.body).getPropertyValue("--tooltip-bg").trim() ||
        (document.body.classList.contains("light") ? "#ffffff" : "#000000");
      const fgColor = getComputedStyle(document.body).getPropertyValue("--tooltip-fg").trim() ||
        (document.body.classList.contains("light") ? "#000000" : "#ffffff");

      let labelText = "Residue " + atom.resi;
      if (currentGene && protData[currentGene] && protData[currentGene].labels) {
        const customLabel = protData[currentGene].labels[atom.resi];
        if (customLabel) labelText = customLabel;
      }

      viewerInstance.addLabel(labelText, {
        position: { x: atom.x, y: atom.y, z: atom.z },
        backgroundColor: bgColor,
        backgroundOpacity: 0.7,
        fontColor: fgColor,
        fontSize: 12,
        showBackground: true,
        inFront: true
      });
    },
    function(atom, viewerInstance) {
      viewerInstance.removeAllLabels();
    }
  );
}

// === Apply Style ===
function applyStyle(gene, style) {
  if (!protData || !protData[gene]) return;

  viewer.clear();
  viewer.addModel(protData[gene].pdb, "pdb");
  enableResidueHoverLabels();

  setTimeout(function() {
    switch (style) {
      case "cartoon":
        viewer.setStyle({}, { cartoon: { color: "spectrum" } });
        break;
      case "surface-plddt":
        viewer.setStyle({}, { cartoon: { color: "lightgrey", opacity: 1.0, thickness: 1.2 } });
        colorSurface(protData[gene].plddt, false, true, "plddt");
        break;
      case "surface-fpocket":
        viewer.setStyle({}, { cartoon: { color: "lightgrey", opacity: 1.0, thickness: 1.2 } });
        colorSurface(protData[gene].fpocket, false, true, "fpocket");
        break;
      case "surface-sasa":
        viewer.setStyle({}, { cartoon: { color: "lightgrey", opacity: 1.0, thickness: 1.2 } });
        colorSurface(protData[gene].sasa, false, true, "sasa");
        break;
      case "surface-pi":
        viewer.setStyle({}, { cartoon: { color: "lightgrey", opacity: 1.0, thickness: 1.2 } });
        colorSurface(protData[gene].pi, true, true, "pi");
        break;
    }

    updateLegend(gene, style);
    viewer.zoomTo();
    viewer.render();
  }, 50);
}

// === DOMContentLoaded ===
document.addEventListener("DOMContentLoaded", function() {
  const cssBg = getComputedStyle(document.body).getPropertyValue("--bg-panel").trim() || "#ffffff";
  viewer = $3Dmol.createViewer("viewer-container", { backgroundColor: cssBg });

  const protSelect = document.getElementById("proteinSelect");
  Object.keys(protData).forEach(gene => {
    const opt = document.createElement("option");
    opt.value = gene;
    opt.textContent = gene;
    protSelect.appendChild(opt);
  });

  if (Object.keys(protData).length > 0) {
    currentGene = Object.keys(protData)[0];
    protSelect.value = currentGene;
    applyStyle(currentGene, "cartoon");
  }

  document.getElementById("styleSelect").addEventListener("change", function(e) {
    if (currentGene) applyStyle(currentGene, e.target.value);
  });

  protSelect.addEventListener("change", function(e) {
    currentGene = e.target.value;
    const style = document.getElementById("styleSelect").value;
    applyStyle(currentGene, style);
  });

  window.addEventListener("message", function(event) {
    if (!event.data || event.data.type !== "theme") return;
    const theme = event.data.theme;
    document.body.classList.toggle("light", theme === "light");
    if (viewer) {
      viewer.setBackgroundColor(theme === "light" ? "#ffffff" : "#0e1430");
      viewer.render();
    }
  });
});
