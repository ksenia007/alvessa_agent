var currentGene = null;
let viewer = null;

// === Color Maps ===
const viridis = [
  [68, 1, 84], [71, 44, 122], [59, 81, 139], [44, 113, 142],
  [33, 144, 141], [39, 173, 129], [92, 200, 99], [170, 220, 50],
  [253, 231, 37]
];
const rdbu = [
  [5, 48, 97], [33, 102, 172], [67, 147, 195], [146, 197, 222],
  [247, 247, 247], [244, 165, 130], [214, 96, 77], [178, 24, 43], [103, 0, 31]
];

// Warm palette for MoRF hotspot:
const hotspotWarm = [
  [180, 220, 180],   // 0.0 soft pale green
  [200, 235, 200],   // 0.25 muted green
  [220, 245, 220],   // 0.5 very pale green (blend with bg)
  [245, 220, 160],   // 0.75 soft orange
  [255, 120, 40]     // 1.0 vivid orange-red hotspot
];

// === Centralized opacity lookup ===
const SURFACE_OPACITY = {
  "plddt": 0.8, "fpocket": 0.8, "sasa": 1.0,
  "pi": 1.0, "disorder": 0.85, "morf": 0.85
};

// === Utility Functions ===
function lerpColor(palette, value) {
  value = Math.min(1, Math.max(0, value || 0));
  const idx = value * (palette.length - 1);
  let i = Math.floor(idx);
  let f = idx - i;
  if (i >= palette.length - 1) { i = palette.length - 1; f = 0; }
  const next = palette[i + 1] || palette[i];
  const rgb = [0, 1, 2].map(j => Math.round(palette[i][j] + f * (next[j] - palette[i][j])));
  return (rgb[0] << 16) | (rgb[1] << 8) | rgb[2];
}
const viridisColor    = v => lerpColor(viridis, v);
const rdbuColor       = v => lerpColor(rdbu, v);
const hotspotWarmColor= v => lerpColor(hotspotWarm, v);

function setHoverLabels(labelFunc) {
  viewer.setHoverable({}, true,
    (atom, v) => {
      const labelText = labelFunc(atom);
      if (labelText) {
        v.addLabel(labelText, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          backgroundColor: "black", backgroundOpacity: 0.7,
          fontColor: "white", fontSize: 12,
          showBackground: true, inFront: true
        });
      }
    },
    (_, v) => v.removeAllLabels()
  );
}

// === PubChem Lookup (Title only) ===
function fetchPubChemName(inchikey, labelElem, valueElem) {
  const block = document.getElementById("ligandNameBlock");
  if (!inchikey) {
    block.style.display = "none";
    return;
  }

  // Show the block immediately with a loading placeholder
  labelElem.style.display = "inline";
  valueElem.textContent = "Loading…";
  valueElem.style.display = "inline";
  block.style.display = "block";

  fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/${inchikey}/cids/JSON`)
    .then(res => res.json())
    .then(data => {
      if (!data.IdentifierList || !data.IdentifierList.CID) {
        block.style.display = "none";
        return;
      }
      const cid = data.IdentifierList.CID[0];
      return fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/Title/JSON`)
        .then(res => res.json())
        .then(propData => {
          const title = propData?.PropertyTable?.Properties?.[0]?.Title;
          if (title) {
            valueElem.textContent = title;
            block.style.display = "block";
          } else {
            block.style.display = "none";
          }
        });
    })
    .catch(err => {
      console.error("PubChem fetch error:", err);
      block.style.display = "none";
    });
}

// === Ligand Dropdown ===
function populateLigandDropdown(gene, ligandSelect, ligandLabel, ligandNameLabel, ligandNameValue) {
  const ligandNameBlock = document.getElementById("ligandNameBlock");
  ligandSelect.innerHTML = "";
  const ligs = protData[gene]?.biolip2;
  if (ligs && ligs.length) {
    ligs.forEach((ligand, idx) => {
      const opt = document.createElement("option");
      opt.value = idx;
      opt.textContent = `Site: ${ligand.bs_code} | Chain: ${ligand.chain_id} | ${ligand.chembl_id} | PDB: ${ligand.pdb_id}`;
      ligandSelect.appendChild(opt);
    });
    ligandSelect.selectedIndex = 0;
    renderLigand(gene, ligs[0]);
    ligandLabel.style.display = "inline";
    ligandSelect.style.display = "inline";
    // Name block visibility handled by fetchPubChemName
    fetchPubChemName(ligs[0].inchikey, ligandNameLabel, ligandNameValue);
  } else {
    ligandLabel.style.display = "none";
    ligandSelect.style.display = "none";
    ligandNameBlock.style.display = "none";
  }
}

// === Legend Drawing ===
function drawLegend(colorFunc) {
  const canvas = document.getElementById("colorLegend");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  for (let i = 0; i < canvas.width; i++) {
    const rgbInt = colorFunc(i / (canvas.width - 1));
    ctx.fillStyle = "#" + rgbInt.toString(16).padStart(6, "0");
    ctx.fillRect(i, 0, 1, canvas.height);
  }
}

function updateLegend(gene, style) {
  const legend = document.getElementById("legend");
  if (!legend) return;
  const minLabel = document.getElementById("legendMin");
  const maxLabel = document.getElementById("legendMax");

  legend.style.display = "none";
  const stats = protData[gene]?.stats;

  const show = (min, max, colorFunc) => {
    minLabel.innerHTML = min; maxLabel.innerHTML = max;
    legend.style.display = "flex"; drawLegend(colorFunc);
  };

  switch (style) {
    case "surface-plddt": show("<strong>Low</strong> (0&#37;)", "<strong>High</strong> (100&#37;)", viridisColor); break;
    case "surface-fpocket": if (stats?.fpocket) show(`<strong>Min</strong> (${stats.fpocket.min.toFixed(3)})`, `<strong>Max</strong> (${stats.fpocket.max.toFixed(3)})`, viridisColor); break;
    case "surface-sasa": if (stats?.sasa) show(`<strong>Min</strong> (${stats.sasa.min.toFixed(2)} A<sup>2</sup>/res)`, `<strong>Max</strong> (${stats.sasa.max.toFixed(2)} A<sup>2</sup>/res)`, viridisColor); break;
    case "surface-pi": if (stats?.pi) show("<strong>-1.0</strong> (apolar)", "<strong>+1.0</strong> (polar)", rdbuColor); break;
    case "surface-disorder": show("<strong>0.0</strong> (ordered)", "<strong>1.0</strong> (disordered)", viridisColor); break;
    case "surface-morf": show("<strong>0.0</strong>", "<strong>1.0</strong> (hotspot)", hotspotWarmColor); break;
  }
}

// === Surface Coloring ===
function colorSurface(contributions, useDiverging = false, mode = "generic") {
  if (!contributions?.length) return;
  const scoreMap = {};
  contributions.forEach(c => {
    const num = Number(c.score);
    if (Number.isFinite(num) && c.residue_no != null) scoreMap[c.residue_no] = Math.min(Math.max(num, 0), 1);
  });
  viewer.addSurface(mode === "sasa" ? $3Dmol.SurfaceType.SAS : $3Dmol.SurfaceType.VDW, {
    transparent: (SURFACE_OPACITY[mode] || 1) < 1.0,
    opacity: SURFACE_OPACITY[mode] || 1.0,
    colorfunc: atom => {
      if (typeof scoreMap[atom?.resi] === "number") {
        if (mode === "pi")   return rdbuColor(scoreMap[atom.resi]);
        if (mode === "morf") return hotspotWarmColor(scoreMap[atom.resi]);
        return viridisColor(scoreMap[atom.resi]);
      }
      return 0xAAAAAA;
    }
  });
}

// === Ligand Rendering (BioLiP2 mode) ===
function renderLigand(gene, ligand) {
  if (!ligand?.pdb_id) return;
  viewer.clear();
  $3Dmol.download(`pdb:${ligand.pdb_id}`, viewer, {}, function() {
    viewer.setStyle({ hetflag: false }, { cartoon: { color: "lightgrey", opacity: 1.0, thickness: 1.2 } });
    viewer.setStyle({ hetflag: true, resn: "HOH" }, {});
    viewer.setStyle({ hetflag: true, not: { resn: "HOH" } }, { wireframe: { radius: 0.01, color: "lightgrey" } });
    const sel = { hetflag: true, resn: ligand.ligand_comp_id, chain: ligand.ligand_chain_id };
    viewer.setStyle(sel, { stick: { radius: 0.25, colorscheme: "Jmol" }, sphere: { scale: 0.35, colorscheme: "Jmol" } });
    viewer.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.5, color: "deepskyblue" }, sel);
    setHoverLabels(atom => {
      if (!atom) return "";
      if (atom.hetflag && atom.resn !== "HOH")
        return `Ligand: ${ligand.ligand_comp_id} | ProtChain: ${ligand.chain_id} | LigChain: ${ligand.ligand_chain_id} | ${ligand.chembl_id} | Atom: ${atom.elem} ${atom.serial}`;
      if (!atom.hetflag) return `ProtChain: ${atom.chain} | ${atom.resn}${atom.resi}`;
      return "";
    });
    viewer.zoomTo(sel); viewer.zoom(0.67, 1000); viewer.render();
  });
}

// === Apply Style ===
function applyStyle(gene, style) {
  if (!protData[gene]) return;
  viewer.clear();
  if (style === "surface-biolip2") {
    updateLegend(gene, "");
    if (protData[gene]?.biolip2?.length) {
      renderLigand(gene, protData[gene].biolip2[0]);
    }
  } else {
    viewer.addModel(protData[gene].pdb, "pdb");
    setHoverLabels(atom => (!atom?.hetflag ? `ProtChain: ${atom.chain} | ${atom.resn}${atom.resi}` : ""));
    setTimeout(() => {
      const setCartoon = () => viewer.setStyle({}, { cartoon: { color: "lightgrey" } });
      switch (style) {
        case "cartoon":        viewer.setStyle({}, { cartoon: { color: "spectrum" } }); break;
        case "surface-plddt":  setCartoon(); colorSurface(protData[gene].plddt,   false, "plddt"); break;
        case "surface-fpocket":setCartoon(); colorSurface(protData[gene].fpocket, false, "fpocket"); break;
        case "surface-sasa":   viewer.setStyle({}, { cartoon: { color: "white", opacity: 0.01 } }); colorSurface(protData[gene].sasa, false, "sasa"); break;
        case "surface-pi":     viewer.setStyle({}, { cartoon: { color: "white", opacity: 0.01 } }); colorSurface(protData[gene].pi,  true,  "pi"); break;
        case "surface-disorder": setCartoon(); colorSurface(protData[gene].disorder, false, "disorder"); break;
        case "surface-morf":     setCartoon(); colorSurface(protData[gene].morf,     false, "morf"); break;
      }
      updateLegend(gene, style); viewer.zoomTo(); viewer.render();
    }, 50);
  }
}

// === DOMContentLoaded ===
document.addEventListener("DOMContentLoaded", () => {
  const cssBg = getComputedStyle(document.body).getPropertyValue("--bg-panel").trim() || "#ffffff";
  viewer = $3Dmol.createViewer("viewer-container", { backgroundColor: cssBg });

  const protSelect = document.getElementById("proteinSelect");
  const styleSelect = document.getElementById("styleSelect");
  const ligandSelect = document.getElementById("ligandSelect");
  const ligandLabel = document.getElementById("ligandLabel");
  const ligandNameLabel = document.getElementById("ligandNameLabel");
  const ligandNameValue = document.getElementById("ligandNameValue");
  const ligandNameBlock = document.getElementById("ligandNameBlock");

  Object.keys(protData).forEach(gene => {
    const opt = document.createElement("option");
    opt.value = gene; opt.textContent = gene;
    protSelect.appendChild(opt);
  });

  if (Object.keys(protData).length) {
    currentGene = Object.keys(protData)[0];
    protSelect.value = currentGene;
    const biolip2Option = Array.from(styleSelect.options).find(opt => opt.value === "surface-biolip2");
    if (biolip2Option) {
      if (protData[currentGene]?.biolip2?.length) {
        biolip2Option.style.display = "block";
      } else {
        biolip2Option.style.display = "none";
        if (styleSelect.value === "surface-biolip2") {
          const firstVisible = Array.from(styleSelect.options).find(o => o.style.display !== "none");
          if (firstVisible) styleSelect.value = firstVisible.value;
        }
      }
    }
    applyStyle(currentGene, styleSelect.value);
  }

  styleSelect.addEventListener("change", e => {
    if (!currentGene) return;
    const style = e.target.value;
    applyStyle(currentGene, style);
    if (style === "surface-biolip2") {
      populateLigandDropdown(currentGene, ligandSelect, ligandLabel, ligandNameLabel, ligandNameValue);
    } else {
      ligandLabel.style.display = "none";
      ligandSelect.style.display = "none";
      ligandNameBlock.style.display = "none";
    }
  });

  ligandSelect.addEventListener("change", e => {
    if (currentGene && protData[currentGene]?.biolip2?.length) {
      const idx = parseInt(e.target.value, 10);
      if (!isNaN(idx)) {
        const ligand = protData[currentGene].biolip2[idx];
        renderLigand(currentGene, ligand);
        fetchPubChemName(ligand.inchikey, ligandNameLabel, ligandNameValue);
      }
    } else {
      ligandNameBlock.style.display = "none";
    }
  });

  protSelect.addEventListener("change", e => {
    currentGene = e.target.value;
    let style = styleSelect.value;
    const biolip2Option = Array.from(styleSelect.options).find(opt => opt.value === "surface-biolip2");
    const hasLigands = !!(protData[currentGene]?.biolip2?.length);
    if (biolip2Option) {
      biolip2Option.style.display = hasLigands ? "block" : "none";
      if (!hasLigands && style === "surface-biolip2") {
        const firstVisible = Array.from(styleSelect.options).find(o => o.style.display !== "none");
        if (firstVisible) {
          styleSelect.value = firstVisible.value;
          style = firstVisible.value;
        }
        ligandLabel.style.display = "none";
        ligandSelect.style.display = "none";
        ligandNameBlock.style.display = "none";
      }
    }
    applyStyle(currentGene, style);
    if (style === "surface-biolip2") {
      populateLigandDropdown(currentGene, ligandSelect, ligandLabel, ligandNameLabel, ligandNameValue);
    }
  });

  window.addEventListener("message", event => {
    if (!event.data || event.data.type !== "theme") return;
    const theme = event.data.theme;
    document.body.classList.toggle("light", theme === "light");
    if (viewer) {
      viewer.setBackgroundColor(theme === "light" ? "#ffffff" : "#0e1430");
      viewer.render();
    }
  });
});
