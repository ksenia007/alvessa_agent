// src/tools/chembl/frontend/tool_chembl.js
// Author: Dmitri Kosenkov
// Updated: 2025-10-06
//
// Full revision:
// - Adds collapsible sections (FDA-approved, Clinical, Bioactivity)
// - Adds "Approval Year" column for FDA-approved drugs
// - Synchronizes light/dark theme with parent UI (like tool_prot.js)
// - Keeps all other functionality identical

let currentGene = null;
let RDKit = null; // RDKit.js instance

// =============================
// Helpers: URLs
// =============================
function chemblUrl(chemblId) {
  return `https://www.ebi.ac.uk/chembl/compound_report_card/${chemblId}/`;
}
function pubchemCidUrl(cid) {
  return `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`;
}

// =============================
// PubChem Lookup (CID + Title)
// =============================
async function fetchPubChemInfo(inchikey) {
  if (!inchikey) return null;
  try {
    const cidResp = await fetch(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/${inchikey}/cids/JSON`
    );
    const cidData = await cidResp.json();
    if (!cidData.IdentifierList?.CID?.length) return null;
    const cid = cidData.IdentifierList.CID[0];

    await new Promise((res) => setTimeout(res, 400)); // avoid rate-limit
    const propResp = await fetch(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/Title/JSON`
    );
    const propData = await propResp.json();
    const title = propData?.PropertyTable?.Properties?.[0]?.Title || "";
    return { cid, title };
  } catch (err) {
    console.error("PubChem lookup error:", err);
    return null;
  }
}

// =============================
// Molecule SVG Renderer
// =============================
function renderMol2D(cell, smiles, inchikey) {
  if (!RDKit || !smiles) {
    cell.textContent = "-";
    return;
  }
  try {
    const mol = RDKit.get_mol(smiles);
    if (!mol) throw new Error("RDKit mol creation failed");
    const svg = mol.get_svg();
    mol.delete();

    const wrapper = document.createElement("div");
    wrapper.className = "mol-thumb";
    wrapper.innerHTML = svg;
    const svgElem = wrapper.querySelector("svg");
    svgElem.style.maxWidth = "90px";
    svgElem.style.cursor = "zoom-in";
    svgElem.title = "Click to enlarge";
    wrapper.addEventListener("click", () => openMolModal(svg, inchikey, smiles));
    cell.innerHTML = "";
    cell.appendChild(wrapper);
  } catch (err) {
    console.warn("RDKit render failed:", err);
    cell.textContent = "-";
  }
}

// =============================
// Modal Viewer
// =============================
async function openMolModal(svg, inchikey, smiles) {
  const modal = document.createElement("div");
  modal.className = "mol-modal";
  modal.innerHTML = `
    <div class="mol-modal-backdrop"></div>
    <div class="mol-modal-content">
      <div class="mol-modal-header">
        <span>Molecule View</span>
        <button class="mol-modal-close">×</button>
      </div>
      <div class="mol-modal-body">${svg}</div>
      <div class="mol-modal-footer">
        <span id="molName" style="flex:1;text-align:left;overflow:hidden;white-space:nowrap;text-overflow:ellipsis;">Loading PubChem…</span>
        <button id="saveSvgBtn">Save SVG</button>
      </div>
    </div>`;
  document.body.appendChild(modal);
  requestAnimationFrame(() => modal.classList.add("active"));
  modal.querySelector(".mol-modal-close").onclick = () => modal.remove();
  modal.querySelector(".mol-modal-backdrop").onclick = () => modal.remove();

  modal.querySelector("#saveSvgBtn").onclick = () => {
    const blob = new Blob([svg], { type: "image/svg+xml" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${inchikey || "molecule"}.svg`;
    a.click();
  };

  const info = await fetchPubChemInfo(inchikey);
  const molNameElem = modal.querySelector("#molName");
  if (info) {
    const { cid, title } = info;
    if (title?.length) {
      let displayTitle = title;
      const MAX_LEN = 60;
      if (displayTitle.length > MAX_LEN)
        displayTitle = displayTitle.slice(0, MAX_LEN - 3) + "…";
      molNameElem.innerHTML = `<a href="${pubchemCidUrl(cid)}" target="_blank" rel="noopener" title="${title}">${displayTitle}</a>`;
    } else if (cid) {
      molNameElem.innerHTML = `<a href="${pubchemCidUrl(cid)}" target="_blank" rel="noopener">PubChem CID ${cid}</a>`;
    } else molNameElem.style.display = "none";
  } else molNameElem.style.display = "none";
}

// =============================
// Table Builder
// =============================
async function buildTable(tableBody, rows, type) {
  tableBody.innerHTML = "";
  if (!rows?.length) return;

  for (const row of rows) {
    const tr = document.createElement("tr");
    let chemblId = "-",
      mtype = "-",
      phaseOrEvidence = "-",
      inchikey = "",
      smiles = "",
      approvalYear = "";

    if (type === "approved_drugs") {
      [chemblId, mtype, inchikey, smiles, approvalYear] = row;
    } else if (type === "clinical_trials") {
      [chemblId, phaseOrEvidence, mtype, inchikey, smiles] = row;
    } else if (type === "bioactivity") {
      [chemblId, mtype, phaseOrEvidence, inchikey, smiles] = row;
    }

    chemblId = chemblId?.trim?.() || "-";
    mtype = mtype?.trim?.() || "-";
    approvalYear = approvalYear?.toString?.().trim?.() || "";
    phaseOrEvidence =
      typeof phaseOrEvidence === "string"
        ? phaseOrEvidence.trim()
        : phaseOrEvidence ?? "";
    inchikey = inchikey?.trim?.() || "";
    smiles = smiles?.trim?.() || "";

    const tdId = document.createElement("td");
    tdId.innerHTML = chemblId.startsWith("CHEMBL")
      ? `<a href="${chemblUrl(chemblId)}" target="_blank" rel="noopener">${chemblId}</a>`
      : "-";
    tr.appendChild(tdId);

    const tdType = document.createElement("td");
    tdType.textContent = mtype || "-";
    tr.appendChild(tdType);

    if (type === "approved_drugs") {
      const tdYear = document.createElement("td");
      tdYear.textContent = approvalYear || "-";
      tr.appendChild(tdYear);
    }

    if (type === "clinical_trials") {
      const tdPhase = document.createElement("td");
      tdPhase.textContent =
        phaseOrEvidence && phaseOrEvidence !== "0"
          ? `Phase ${phaseOrEvidence}`
          : "Preclinical";
      tr.appendChild(tdPhase);
    } else if (type === "bioactivity") {
      const tdEv = document.createElement("td");
      tdEv.textContent = phaseOrEvidence || "-";
      tr.appendChild(tdEv);
    }

    const tdMol = document.createElement("td");
    smiles ? renderMol2D(tdMol, smiles, inchikey) : (tdMol.textContent = "-");
    tr.appendChild(tdMol);

    tableBody.appendChild(tr);
  }
}

// =============================
// Populate Viewer
// =============================
async function populateGene(gene) {
  if (!chemblData[gene]) return;
  currentGene = gene;
  /*document.getElementById("geneTitle").textContent =
    `ChEMBL Drug-Target Data for ${gene}`;*/

  const sections = [
    { key: "approved_drugs", id: "approvedSection", tbl: "approvedTable", title: "FDA-approved Drugs" },
    { key: "clinical_trials", id: "clinicalSection", tbl: "clinicalTable", title: "Clinical & Preclinical Trials" },
    { key: "bioactivity", id: "bioactivitySection", tbl: "bioactivityTable", title: "Bioactivity Evidence" },
  ];

  for (const { key, id, tbl, title } of sections) {
    const data = chemblData[gene][key];
    const section = document.getElementById(id);
    const content = section.querySelector(".collapsible-content");
    const header = section.querySelector(".collapsible-header .title");
    const count = section.querySelector(".collapsible-header .count");

    if (!data?.length) {
      section.style.display = "none";
      continue;
    }

    header.textContent = title;
    count.textContent = data.length;
    section.style.display = "block";
    await buildTable(document.getElementById(tbl), data, key);
  }
}

// =============================
// DOM Ready + Theme Sync
// =============================
document.addEventListener("DOMContentLoaded", async () => {
  RDKit = await initRDKitModule();
  console.log("RDKit.js ready:", !!RDKit);

  const geneSelect = document.getElementById("geneSelect");
  Object.keys(chemblData).forEach((g) => {
    const opt = document.createElement("option");
    opt.value = g;
    opt.textContent = g;
    geneSelect.appendChild(opt);
  });

  if (Object.keys(chemblData).length) {
    const firstGene = Object.keys(chemblData)[0];
    geneSelect.value = firstGene;
    populateGene(firstGene);
  }

  // Collapsible section toggles
  document.querySelectorAll(".collapsible-header").forEach((header) => {
    const section = header.closest(".collapsible-section");
    const content = section.querySelector(".collapsible-content");
    const icon = header.querySelector(".toggle-icon");

    header.addEventListener("click", () => {
      const isOpen = content.classList.toggle("open");
      header.classList.toggle("active", isOpen);
      icon.textContent = isOpen ? "▾" : "▸";
    });

    // start collapsed
    content.classList.remove("open");
    icon.textContent = "▸";
  });

  geneSelect.addEventListener("change", (e) => populateGene(e.target.value));

  // --- Theme Sync (identical to tool_prot.js) ---
  try {
    const parentIsLight =
      window.parent?.document?.body?.classList.contains("light");
    if (parentIsLight) document.body.classList.add("light");
  } catch (err) {
    console.warn("Theme sync: parent detection failed", err);
  }

  window.addEventListener("message", (event) => {
    if (!event.data || event.data.type !== "theme") return;
    const theme = event.data.theme;
    document.body.classList.toggle("light", theme === "light");

    // Smooth transition for background and text colors
    const els = document.querySelectorAll(
      ".collapsible-header, .chembl-table, body"
    );
    els.forEach((el) => {
      el.style.transition = "background-color 0.25s ease, color 0.25s ease";
      if (theme === "light") {
        el.style.background = "#ffffff";
        el.style.color = "#0b1020";
      } else {
        el.style.background = "#1b2238";
        el.style.color = "#f2f5ff";
      }
    });
  });
});
