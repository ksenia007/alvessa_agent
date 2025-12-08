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
  [180, 220, 180],
  [200, 235, 200],
  [220, 245, 220],
  [245, 220, 160],
  [255, 120, 40]
];

// CysDB overlay color (high contrast vs viridis)
const CYSDB_COLOR = 0xff1f8f;

// === Centralized opacity lookup ===
const SURFACE_OPACITY = {
  "plddt": 0.8,
  "fpocket": 0.8,
  "sasa": 1.0,
  "pi": 1.0,
  "disorder": 0.85,
  "morf": 0.85
};

function lerpColor(palette, value) {
  value = Math.min(1, Math.max(0, value || 0));
  const idx = value * (palette.length - 1);
  let i = Math.floor(idx);
  let f = idx - i;
  if (i >= palette.length - 1) {
    i = palette.length - 1;
    f = 0;
  }
  const next = palette[i + 1] || palette[i];
  const rgb = [0, 1, 2].map(function(j) {
    return Math.round(palette[i][j] + f * (next[j] - palette[i][j]));
  });
  return (rgb[0] << 16) | (rgb[1] << 8) | rgb[2];
}
const viridisColor = function(v) { return lerpColor(viridis, v); };
const rdbuColor = function(v) { return lerpColor(rdbu, v); };
const hotspotWarmColor = function(v) { return lerpColor(hotspotWarm, v); };

function setHoverLabels(labelFunc) {
  viewer.setHoverable(
    {},
    true,
    function(atom, v) {
      const labelText = labelFunc(atom);
      if (labelText) {
        v.addLabel(labelText, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          backgroundColor: "black",
          backgroundOpacity: 0.7,
          fontColor: "white",
          fontSize: 12,
          showBackground: true,
          inFront: true
        });
      }
    },
    function(_, v) { v.removeAllLabels(); }
  );
}

// === CysDB helpers ===
function getCysdbTracks(gene) {
  const data = protData[gene] || {};
  return {
    hyperreactive: data.cysdb_hyperreactive || [],
    ligandable: data.cysdb_ligandable || [],
    isActSite: data.cysdb_is_act_site || [],
    nearActSite: data.cysdb_near_act_site || [],
    isBindSite: data.cysdb_is_bind_site || [],
    nearBindSite: data.cysdb_near_bind_site || []
  };
}

function hasAnyFlag(track) {
  if (!track || !track.length) return false;
  for (let i = 0; i < track.length; i++) {
    if (Number(track[i].score) === 1) return true;
  }
  return false;
}

function geneHasAnyCysdb(gene) {
  if (!gene || !protData[gene]) return false;
  const tracks = getCysdbTracks(gene);
  return (
    hasAnyFlag(tracks.hyperreactive) ||
    hasAnyFlag(tracks.ligandable) ||
    hasAnyFlag(tracks.isActSite) ||
    hasAnyFlag(tracks.nearActSite) ||
    hasAnyFlag(tracks.isBindSite) ||
    hasAnyFlag(tracks.nearBindSite)
  );
}

function updateFpocketLabel(gene) {
  const styleSelect = document.getElementById("styleSelect");
  if (!styleSelect || !gene) return;
  const opt = Array.prototype.slice.call(styleSelect.options).find(function(o) {
    return o.value === "surface-fpocket";
  });
  if (!opt) return;
  if (geneHasAnyCysdb(gene)) {
    opt.textContent = "Druggablity (FPocket+CysDB)";
  } else {
    opt.textContent = "Druggablity (FPocket)";
  }
}

function updateCysdbControls(gene, style) {
  const cysdbControls = document.getElementById("cysdbControls");
  if (!cysdbControls) return;

  const fpocketCb = document.getElementById("cysdbFpocketBase");

  // Only relevant in FPocket mode and when this gene exists
  if (!gene || style !== "surface-fpocket" || !protData[gene]) {
    if (fpocketCb) {
      fpocketCb.checked = true;
      fpocketCb.disabled = true;
    }
    cysdbControls.style.display = "none";
    return;
  }

  const tracks = getCysdbTracks(gene);
  const hasHyper = hasAnyFlag(tracks.hyperreactive);
  const hasLigandable = hasAnyFlag(tracks.ligandable);
  const hasIsAct = hasAnyFlag(tracks.isActSite);
  const hasNearAct = hasAnyFlag(tracks.nearActSite);
  const hasIsBind = hasAnyFlag(tracks.isBindSite);
  const hasNearBind = hasAnyFlag(tracks.nearBindSite);

  const hasIdentified = hasHyper || hasLigandable || hasIsAct || hasNearAct || hasIsBind || hasNearBind;
  const anyFlag = geneHasAnyCysdb(gene);

  // If no CysDB flags at all for this gene, hide the entire line
  if (!anyFlag) {
    if (fpocketCb) {
      fpocketCb.checked = true;
      fpocketCb.disabled = true;
    }
    cysdbControls.style.display = "none";
    return;
  }

  // Show line and enable FPocket base checkbox
  cysdbControls.style.display = "flex";
  if (fpocketCb) {
    fpocketCb.disabled = false;
    fpocketCb.checked = true;
  }

  function setCheckbox(id, enabled) {
    const cb = document.getElementById(id);
    if (!cb) return;
    cb.disabled = !enabled;
    cb.checked = enabled;
  }

  setCheckbox("cysdbIdentified", hasIdentified);
  setCheckbox("cysdbHyperreactive", hasHyper);
  setCheckbox("cysdbLigandable", hasLigandable);
  setCheckbox("cysdbIsActSite", hasIsAct);
  setCheckbox("cysdbNearActSite", hasNearAct);
  setCheckbox("cysdbIsBindSite", hasIsBind);
  setCheckbox("cysdbNearBindSite", hasNearBind);
}

function buildCysdbMaskForCurrentState(gene) {
  if (!gene || !protData[gene]) return null;

  const cysdbControls = document.getElementById("cysdbControls");
  if (!cysdbControls || cysdbControls.style.display === "none") return null;

  const tracks = getCysdbTracks(gene);
  const mask = {};

  const cbIdentified = document.getElementById("cysdbIdentified");
  const cbHyper = document.getElementById("cysdbHyperreactive");
  const cbLigandable = document.getElementById("cysdbLigandable");
  const cbIsAct = document.getElementById("cysdbIsActSite");
  const cbNearAct = document.getElementById("cysdbNearActSite");
  const cbIsBind = document.getElementById("cysdbIsBindSite");
  const cbNearBind = document.getElementById("cysdbNearBindSite");

  function applyTrack(track) {
    if (!track || !track.length) return;
    for (let i = 0; i < track.length; i++) {
      const entry = track[i];
      if (Number(entry.score) === 1 && entry.residue_no != null) {
        mask[entry.residue_no] = true;
      }
    }
  }

  const identifiedOn = cbIdentified && !cbIdentified.disabled && cbIdentified.checked;

  if (identifiedOn) {
    applyTrack(tracks.hyperreactive);
    applyTrack(tracks.ligandable);
    applyTrack(tracks.isActSite);
    applyTrack(tracks.nearActSite);
    applyTrack(tracks.isBindSite);
    applyTrack(tracks.nearBindSite);
  } else {
    if (cbHyper && !cbHyper.disabled && cbHyper.checked) {
      applyTrack(tracks.hyperreactive);
    }
    if (cbLigandable && !cbLigandable.disabled && cbLigandable.checked) {
      applyTrack(tracks.ligandable);
    }
    if (cbIsAct && !cbIsAct.disabled && cbIsAct.checked) {
      applyTrack(tracks.isActSite);
    }
    if (cbNearAct && !cbNearAct.disabled && cbNearAct.checked) {
      applyTrack(tracks.nearActSite);
    }
    if (cbIsBind && !cbIsBind.disabled && cbIsBind.checked) {
      applyTrack(tracks.isBindSite);
    }
    if (cbNearBind && !cbNearBind.disabled && cbNearBind.checked) {
      applyTrack(tracks.nearBindSite);
    }
  }

  if (!Object.keys(mask).length) return null;
  return mask;
}

// === PubChem Lookup (Title only) ===
function fetchPubChemName(inchikey, labelElem, valueElem) {
  const block = document.getElementById("ligandNameBlock");
  if (!inchikey) {
    block.style.display = "none";
    return;
  }

  labelElem.style.display = "inline";
  valueElem.textContent = "Loading…";
  valueElem.style.display = "inline";
  block.style.display = "block";

  fetch("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/" + inchikey + "/cids/JSON")
    .then(function(res) { return res.json(); })
    .then(function(data) {
      if (!data.IdentifierList || !data.IdentifierList.CID) {
        block.style.display = "none";
        return;
      }
      const cid = data.IdentifierList.CID[0];
      return fetch("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + cid + "/property/Title/JSON")
        .then(function(res) { return res.json(); })
        .then(function(propData) {
          const title =
            propData &&
            propData.PropertyTable &&
            propData.PropertyTable.Properties &&
            propData.PropertyTable.Properties[0] &&
            propData.PropertyTable.Properties[0].Title;
          if (title) {
            valueElem.textContent = title;
            block.style.display = "block";
          } else {
            block.style.display = "none";
          }
        });
    })
    .catch(function(err) {
      console.error("PubChem fetch error:", err);
      block.style.display = "none";
    });
}

// === Ligand Dropdown ===
function populateLigandDropdown(gene, ligandSelect, ligandLabel, ligandNameLabel, ligandNameValue) {
  const ligandNameBlock = document.getElementById("ligandNameBlock");
  ligandSelect.innerHTML = "";
  const ligs = protData[gene] && protData[gene].biolip2;
  if (ligs && ligs.length) {
    ligs.forEach(function(ligand, idx) {
      const opt = document.createElement("option");
      opt.value = idx;
      opt.textContent =
        "Site: " +
        ligand.bs_code +
        " | Chain: " +
        ligand.chain_id +
        " | " +
        ligand.chembl_id +
        " | PDB: " +
        ligand.pdb_id;
      ligandSelect.appendChild(opt);
    });
    ligandSelect.selectedIndex = 0;
    renderLigand(gene, ligs[0]);
    ligandLabel.style.display = "inline";
    ligandSelect.style.display = "inline";
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
  const stats = protData[gene] && protData[gene].stats;

  function show(min, max, colorFunc) {
    minLabel.innerHTML = min;
    maxLabel.innerHTML = max;
    legend.style.display = "flex";
    drawLegend(colorFunc);
  }

  switch (style) {
    case "surface-plddt":
      show(
        "<strong>Low</strong> (0&#37;)",
        "<strong>High</strong> (100&#37;)",
        viridisColor
      );
      break;
    case "surface-fpocket":
      if (stats && stats.fpocket) {
        show(
          "<strong>Min</strong> (" + stats.fpocket.min.toFixed(3) + ")",
          "<strong>Max</strong> (" + stats.fpocket.max.toFixed(3) + ")",
          viridisColor
        );
      }
      break;
    case "surface-sasa":
      if (stats && stats.sasa) {
        show(
          "<strong>Min</strong> (" + stats.sasa.min.toFixed(2) + " A<sup>2</sup>/res)",
          "<strong>Max</strong> (" + stats.sasa.max.toFixed(2) + " A<sup>2</sup>/res)",
          viridisColor
        );
      }
      break;
    case "surface-pi":
      if (stats && stats.pi) {
        show(
          "<strong>-1.0</strong> (apolar)",
          "<strong>+1.0</strong> (polar)",
          rdbuColor
        );
      }
      break;
    case "surface-disorder":
      show(
        "<strong>0.0</strong> (ordered)",
        "<strong>1.0</strong> (disordered)",
        viridisColor
      );
      break;
    case "surface-morf":
      show("<strong>0.0</strong>", "<strong>1.0</strong> (hotspot)", hotspotWarmColor);
      break;
  }
}

// === Surface Coloring ===
function colorSurface(contributions, useDiverging, mode, cysdbMask, baseEnabled) {
  if (!contributions || !contributions.length) return;
  if (useDiverging === undefined) useDiverging = false;
  if (mode === undefined) mode = "generic";

  // NEW: if FPocket checkbox is off, do not draw FPocket surface at all
  if (mode === "fpocket" && baseEnabled === false) {
    return;
  }

  const scoreMap = {};
  contributions.forEach(function(c) {
    const num = Number(c.score);
    if (Number.isFinite(num) && c.residue_no != null) {
      scoreMap[c.residue_no] = Math.min(Math.max(num, 0), 1);
    }
  });

  viewer.addSurface(
    mode === "sasa" ? $3Dmol.SurfaceType.SAS : $3Dmol.SurfaceType.VDW,
    {
      transparent: (SURFACE_OPACITY[mode] || 1) < 1.0,
      opacity: SURFACE_OPACITY[mode] || 1.0,
      colorfunc: function(atom) {
        if (!atom) return 0xaaaaaa;

        const resi = atom.resi;

        // For FPocket mode, CysDB mask overrides surface color
        if (mode === "fpocket" && cysdbMask && cysdbMask[resi]) {
          return CYSDB_COLOR;
        }

        const score = scoreMap[resi];
        if (typeof score === "number") {
          if (mode === "pi") {
            return rdbuColor(score);
          }
          if (mode === "morf") {
            return hotspotWarmColor(score);
          }
          return viridisColor(score);
        }
        return 0xaaaaaa;
      }
    }
  );
}

// === Ligand Rendering (BioLiP2 mode) ===
function renderLigand(gene, ligand) {
  if (!ligand || !ligand.pdb_id) return;
  viewer.clear();
  $3Dmol.download("pdb:" + ligand.pdb_id, viewer, {}, function() {
    viewer.setStyle(
      { hetflag: false },
      { cartoon: { color: "lightgrey", opacity: 1.0, thickness: 1.2 } }
    );
    viewer.setStyle({ hetflag: true, resn: "HOH" }, {});
    viewer.setStyle(
      { hetflag: true, not: { resn: "HOH" } },
      { wireframe: { radius: 0.01, color: "lightgrey" } }
    );
    const sel = {
      hetflag: true,
      resn: ligand.ligand_comp_id,
      chain: ligand.ligand_chain_id
    };
    viewer.setStyle(sel, {
      stick: { radius: 0.25, colorscheme: "Jmol" },
      sphere: { scale: 0.35, colorscheme: "Jmol" }
    });
    viewer.addSurface(
      $3Dmol.SurfaceType.VDW,
      { opacity: 0.5, color: "deepskyblue" },
      sel
    );
    setHoverLabels(function(atom) {
      if (!atom) return "";
      if (atom.hetflag && atom.resn !== "HOH") {
        return (
          "Ligand: " +
          ligand.ligand_comp_id +
          " | ProtChain: " +
          ligand.chain_id +
          " | LigChain: " +
          ligand.ligand_chain_id +
          " | " +
          ligand.chembl_id +
          " | Atom: " +
          atom.elem +
          " " +
          atom.serial
        );
      }
      if (!atom.hetflag) {
        return "ProtChain: " + atom.chain + " | " + atom.resn + atom.resi;
      }
      return "";
    });
    viewer.zoomTo(sel);
    viewer.zoom(0.67, 1000);
    viewer.render();
  });
}

// === Apply Style ===
function applyStyle(gene, style) {
  if (!protData[gene]) return;
  viewer.clear();
  if (style === "surface-biolip2") {
    updateLegend(gene, "");
    if (protData[gene].biolip2 && protData[gene].biolip2.length) {
      renderLigand(gene, protData[gene].biolip2[0]);
    }
  } else {
    viewer.addModel(protData[gene].pdb, "pdb");
    setHoverLabels(function(atom) {
      if (!atom || atom.hetflag) return "";
      return "ProtChain: " + atom.chain + " | " + atom.resn + atom.resi;
    });
    setTimeout(function() {
      const setCartoon = function() {
        viewer.setStyle({}, { cartoon: { color: "lightgrey" } });
      };
      switch (style) {
        case "cartoon":
          viewer.setStyle({}, { cartoon: { color: "spectrum" } });
          break;
        case "surface-plddt":
          setCartoon();
          colorSurface(protData[gene].plddt, false, "plddt");
          break;
        case "surface-fpocket": {
          const cysdbMask = buildCysdbMaskForCurrentState(gene);
          const baseCb = document.getElementById("cysdbFpocketBase");
          let fpocketEnabled = true;
          if (baseCb && !baseCb.disabled && !baseCb.checked) {
            fpocketEnabled = false;
          }

          // Cartoon backbone: CysDB overrides grey when mask present
          if (cysdbMask) {
            viewer.setStyle({}, {
              cartoon: {
                colorfunc: function(atom) {
                  if (!atom || atom.hetflag) return 0xd3d3d3;
                  const resi = atom.resi;
                  if (cysdbMask && cysdbMask[resi]) {
                    return CYSDB_COLOR;
                  }
                  return 0xd3d3d3;
                }
              }
            });
          } else {
            viewer.setStyle({}, { cartoon: { color: "lightgrey" } });
          }

          // FPocket surface: suppressed when fpocketEnabled is false
          colorSurface(
            protData[gene].fpocket,
            false,
            "fpocket",
            cysdbMask,
            fpocketEnabled
          );
          break;
        }
        case "surface-sasa":
          viewer.setStyle({}, { cartoon: { color: "white", opacity: 0.01 } });
          colorSurface(protData[gene].sasa, false, "sasa");
          break;
        case "surface-pi":
          viewer.setStyle({}, { cartoon: { color: "white", opacity: 0.01 } });
          colorSurface(protData[gene].pi, true, "pi");
          break;
        case "surface-disorder":
          setCartoon();
          colorSurface(protData[gene].disorder, false, "disorder");
          break;
        case "surface-morf":
          setCartoon();
          colorSurface(protData[gene].morf, false, "morf");
          break;
      }
      updateLegend(gene, style);
      viewer.zoomTo();
      viewer.render();
    }, 50);
  }
}

// === DOMContentLoaded ===
document.addEventListener("DOMContentLoaded", function() {
  const cssBg =
    getComputedStyle(document.body).getPropertyValue("--bg-panel").trim() ||
    "#ffffff";
  viewer = $3Dmol.createViewer("viewer-container", { backgroundColor: cssBg });

  const protSelect = document.getElementById("proteinSelect");
  const styleSelect = document.getElementById("styleSelect");
  const ligandSelect = document.getElementById("ligandSelect");
  const ligandLabel = document.getElementById("ligandLabel");
  const ligandNameLabel = document.getElementById("ligandNameLabel");
  const ligandNameValue = document.getElementById("ligandNameValue");
  const ligandNameBlock = document.getElementById("ligandNameBlock");

  const cysdbCheckboxIds = [
    "cysdbFpocketBase",
    "cysdbIdentified",
    "cysdbHyperreactive",
    "cysdbLigandable",
    "cysdbIsActSite",
    "cysdbNearActSite",
    "cysdbIsBindSite",
    "cysdbNearBindSite"
  ];

  Object.keys(protData).forEach(function(gene) {
    const opt = document.createElement("option");
    opt.value = gene;
    opt.textContent = gene;
    protSelect.appendChild(opt);
  });

  if (Object.keys(protData).length) {
    currentGene = Object.keys(protData)[0];
    protSelect.value = currentGene;
    let style = styleSelect.value;

    const biolip2Option = Array.prototype.slice.call(styleSelect.options).find(
      function(opt) { return opt.value === "surface-biolip2"; }
    );
    if (biolip2Option) {
      if (
        protData[currentGene] &&
        protData[currentGene].biolip2 &&
        protData[currentGene].biolip2.length
      ) {
        biolip2Option.style.display = "block";
      } else {
        biolip2Option.style.display = "none";
        if (style === "surface-biolip2") {
          const firstVisible = Array.prototype.slice
            .call(styleSelect.options)
            .find(function(o) { return o.style.display !== "none"; });
          if (firstVisible) {
            styleSelect.value = firstVisible.value;
            style = firstVisible.value;
          }
        }
      }
    }
    applyStyle(currentGene, style);
    updateCysdbControls(currentGene, style);
    updateFpocketLabel(currentGene);
  }

  styleSelect.addEventListener("change", function(e) {
    if (!currentGene) return;
    const style = e.target.value;
    applyStyle(currentGene, style);
    updateCysdbControls(currentGene, style);

    if (style === "surface-biolip2") {
      populateLigandDropdown(
        currentGene,
        ligandSelect,
        ligandLabel,
        ligandNameLabel,
        ligandNameValue
      );
    } else {
      ligandLabel.style.display = "none";
      ligandSelect.style.display = "none";
      ligandNameBlock.style.display = "none";
    }
  });

  ligandSelect.addEventListener("change", function(e) {
    if (
      currentGene &&
      protData[currentGene] &&
      protData[currentGene].biolip2 &&
      protData[currentGene].biolip2.length
    ) {
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

  protSelect.addEventListener("change", function(e) {
    currentGene = e.target.value;
    let style = styleSelect.value;
    const biolip2Option = Array.prototype.slice.call(styleSelect.options).find(
      function(opt) { return opt.value === "surface-biolip2"; }
    );
    const hasLigands =
      !!(
        protData[currentGene] &&
        protData[currentGene].biolip2 &&
        protData[currentGene].biolip2.length
      );

    if (biolip2Option) {
      biolip2Option.style.display = hasLigands ? "block" : "none";
      if (!hasLigands && style === "surface-biolip2") {
        const firstVisible = Array.prototype.slice
          .call(styleSelect.options)
          .find(function(o) { return o.style.display !== "none"; });
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
    updateCysdbControls(currentGene, style);
    updateFpocketLabel(currentGene);

    if (style === "surface-biolip2") {
      populateLigandDropdown(
        currentGene,
        ligandSelect,
        ligandLabel,
        ligandNameLabel,
        ligandNameValue
      );
    }
  });

  // Recolor FPocket surface when CysDB or FPocket checkboxes change
  cysdbCheckboxIds.forEach(function(id) {
    const cb = document.getElementById(id);
    if (!cb) return;
    cb.addEventListener("change", function() {
      if (!currentGene) return;
      const style = styleSelect.value;
      if (style === "surface-fpocket") {
        applyStyle(currentGene, style);
      }
    });
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
