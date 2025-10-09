// src/tools/chembl/frontend/tool_chembl.js
// Author: Dmitri Kosenkov
// Updated: 2025-10-09 (secure postMessage theme sync, modal-only PubChem lookup)
//
// Columns (approved):
//   ChEMBL ID (pref name + withdrawn under ID) | Molecule Type | First Approval | Indications | Mechanism | 2D
//
// Changes:
// - Group approved rows by (chembl_id, first_approval); combine indications.
// - Show (withdrawn) under ChEMBL ID; DO NOT show black box flag there.
// - Always render a "Black Box Warning" panel if black_box===true, even without text/link.
// - Warning row is shifted right (starts under Molecule Type).
// - Collapsible icons: ">" (closed) / "v" (open).
// - Proper grouped count; HTML-escape; RDKit.js render.
// - Secure theme sync via postMessage (works without allow-same-origin).
// - Visual withdrawn cue via CSS hooks: <tr class="is-withdrawn"> and <a class="id-link">.
// - PubChem lookup only inside modal (interactive), not during table build.

let currentGene = null;
let RDKit = null;

// =============================
// Helpers
// =============================
function chemblUrl(chemblId) {
  return "https://www.ebi.ac.uk/chembl/compound_report_card/" + chemblId + "/";
}
function pubchemCidUrl(cid) {
  return "https://pubchem.ncbi.nlm.nih.gov/compound/" + cid;
}
function esc(s) {
  if (s == null) return "";
  return String(s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;");
}

// =============================
// PubChem Lookup (cached + throttle)
// =============================
const _pubchemCache = new Map(); // inchikey -> { cid, title }

async function fetchPubChemInfo(inchikey) {
  if (!inchikey) return null;
  if (_pubchemCache.has(inchikey)) return _pubchemCache.get(inchikey);
  try {
    await new Promise((res) => setTimeout(res, 200));
    const cidResp = await fetch(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/" +
        inchikey +
        "/cids/JSON"
    );
    if (!cidResp.ok) {
      console.warn("PubChem returned " + cidResp.status + " for " + inchikey);
      _pubchemCache.set(inchikey, null);
      return null;
    }
    const cidData = await cidResp.json();
    if (!cidData.IdentifierList || !cidData.IdentifierList.CID?.length) {
      _pubchemCache.set(inchikey, null);
      return null;
    }
    const cid = cidData.IdentifierList.CID[0];
    await new Promise((res) => setTimeout(res, 200));
    const propResp = await fetch(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" +
        cid +
        "/property/Title/JSON"
    );
    if (!propResp.ok) {
      console.warn("PubChem property returned " + propResp.status + " for " + inchikey);
      _pubchemCache.set(inchikey, null);
      return null;
    }
    const propData = await propResp.json();
    const title = propData?.PropertyTable?.Properties?.[0]?.Title || "";
    const out = { cid, title };
    _pubchemCache.set(inchikey, out);
    return out;
  } catch (err) {
    console.warn("PubChem lookup error:", err, "for", inchikey);
    _pubchemCache.set(inchikey, null);
    return null;
  }
}

// =============================
// RDKit SVG rendering
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
    if (svgElem) {
      svgElem.style.maxWidth = "90px";
      svgElem.style.cursor = "zoom-in";
    }
    wrapper.addEventListener("click", () => openMolModal(svg, inchikey, smiles));
    cell.innerHTML = "";
    cell.appendChild(wrapper);
  } catch (err) {
    console.warn("RDKit render failed:", err);
    cell.textContent = "-";
  }
}

// =============================
// Modal Viewer (with PubChem lookup)
// =============================
async function openMolModal(svg, inchikey, smiles) {
  const modal = document.createElement("div");
  modal.className = "mol-modal";
  modal.innerHTML = `
    <div class="mol-modal-backdrop"></div>
    <div class="mol-modal-content">
      <div class="mol-modal-header">
        <span>Molecule View</span>
        <button class="mol-modal-close" aria-label="Close">X</button>
      </div>
      <div class="mol-modal-body">${svg}</div>
      <div class="mol-modal-footer">
        <span id="molName">Loading PubChem...</span>
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
    a.download = (inchikey || "molecule") + ".svg";
    a.click();
  };

  const info = await fetchPubChemInfo(inchikey);
  const molNameElem = modal.querySelector("#molName");
  if (info) {
    const { cid, title } = info;
    if (title) {
      molNameElem.innerHTML =
        '<a href="' +
        pubchemCidUrl(cid) +
        '" target="_blank" rel="noopener">' +
        esc(title) +
        "</a>";
    } else if (cid) {
      molNameElem.innerHTML =
        '<a href="' +
        pubchemCidUrl(cid) +
        '" target="_blank" rel="noopener">PubChem CID ' +
        cid +
        "</a>";
    } else {
      molNameElem.style.display = "none";
    }
  } else {
    molNameElem.style.display = "none";
  }
}

// =============================
// Grouping helpers (Approved)
// =============================
function groupApprovedRows(rows) {
  const groups = new Map();
  for (const d of rows) {
    const key = d.chembl_id + "||" + (d.first_approval || "");
    if (!groups.has(key)) {
      groups.set(key, {
        chembl_id: d.chembl_id,
        first_approval: d.first_approval || "",
        molecule_type: d.molecule_type || "",
        inchi_key: d.inchi_key || "",
        smiles: d.smiles || "",
        withdrawn: !!d.withdrawn,
        black_box: !!d.black_box,
        black_box_text: d.black_box_text || "",
        black_box_url: d.black_box_url || "",
        pref_name: d.pref_name || "",
        indications: new Set(),
        mechanisms: new Set(),
      });
    }
    const g = groups.get(key);
    if (d.indication) g.indications.add(d.indication);
    if (d.mechanism) g.mechanisms.add(d.mechanism);
    if (d.withdrawn) g.withdrawn = true;
    if (d.black_box) g.black_box = true;
    if (!g.black_box_text && d.black_box_text) g.black_box_text = d.black_box_text;
    if (!g.black_box_url && d.black_box_url) g.black_box_url = d.black_box_url;
    if (!g.smiles && d.smiles) g.smiles = d.smiles;
    if (!g.inchi_key && d.inchi_key) g.inchi_key = d.inchi_key;
    if (!g.pref_name && d.pref_name) g.pref_name = d.pref_name;
  }
  const out = [];
  for (const g of groups.values()) {
    out.push({
      chembl_id: g.chembl_id,
      first_approval: g.first_approval,
      molecule_type: g.molecule_type,
      inchi_key: g.inchi_key,
      smiles: g.smiles,
      withdrawn: g.withdrawn,
      black_box: g.black_box,
      black_box_text: g.black_box_text,
      black_box_url: g.black_box_url,
      pref_name: g.pref_name,
      indication: Array.from(g.indications).sort().join(", "),
      mechanism: Array.from(g.mechanisms).sort().join("; "),
    });
  }
  return out;
}

// =============================
// Table Builders
// =============================
async function buildApprovedTable(tbody, groupedRows) {
  tbody.innerHTML = "";
  if (!groupedRows?.length) return;
  for (const drug of groupedRows) {
    const tr = document.createElement("tr");
    if (drug.withdrawn) tr.classList.add("is-withdrawn");

    const tdId = document.createElement("td");
    const idLink =
      '<a class="id-link" href="' +
      chemblUrl(drug.chembl_id) +
      '" target="_blank" rel="noopener">' +
      esc(drug.chembl_id) +
      "</a>";
    const flags = drug.withdrawn ? "(withdrawn)" : "";
    tdId.innerHTML =
      "<div>" +
      idLink +
      "</div>" +
      (drug.pref_name
        ? '<div class="sub muted">' + esc(drug.pref_name) + "</div>"
        : "") +
      (flags ? '<div class="sub withdrawn">' + esc(flags) + "</div>" : "");
    tr.appendChild(tdId);

    const tdType = document.createElement("td");
    tdType.textContent = drug.molecule_type || "-";
    tr.appendChild(tdType);

    const tdYear = document.createElement("td");
    tdYear.textContent = drug.first_approval || "-";
    tr.appendChild(tdYear);

    const tdInd = document.createElement("td");
    tdInd.textContent = drug.indication || "-";
    tr.appendChild(tdInd);

    const tdMoA = document.createElement("td");
    tdMoA.textContent = drug.mechanism || "-";
    tr.appendChild(tdMoA);

    const tdMol = document.createElement("td");
    drug.smiles ? renderMol2D(tdMol, drug.smiles, drug.inchi_key) : (tdMol.textContent = "-");
    tr.appendChild(tdMol);

    tbody.appendChild(tr);

    if (drug.black_box || drug.black_box_text) {
      const warnTr = document.createElement("tr");
      warnTr.className = "bb-warning-row";
      warnTr.innerHTML =
        '<td></td>' +
        '<td colspan="4">' +
        '<div class="bb-warning-box">' +
        '<div class="bb-warning-header">Black Box Warning</div>' +
        '<div class="bb-warning-text">' +
        (esc(drug.black_box_text) ||
          "<em>(No summary text was found from openFDA for this label.)</em>") +
        "</div>" +
        (drug.black_box_url
          ? '<div class="bb-warning-link"><a href="' +
            esc(drug.black_box_url) +
            '" target="_blank" rel="noopener">View FDA Label</a></div>'
          : "") +
        "</div>" +
        "</td>" +
        "<td></td>";
      tbody.appendChild(warnTr);
    }
  }
}

async function buildTable(tbody, rows, type) {
  tbody.innerHTML = "";
  if (!rows?.length) return;
  for (const row of rows) {
    const tr = document.createElement("tr");
    if (type === "clinical_trials") {
      const [chemblId, phase, mtype, inchikey, smiles] = row;
      const tdId = document.createElement("td");
      tdId.innerHTML =
        '<a href="' +
        chemblUrl(chemblId) +
        '" target="_blank" rel="noopener">' +
        esc(chemblId) +
        "</a>";
      tr.appendChild(tdId);
      const tdType = document.createElement("td");
      tdType.textContent = mtype || "-";
      tr.appendChild(tdType);
      const tdPhase = document.createElement("td");
      tdPhase.textContent = phase === 0 ? "Preclinical" : "Phase " + phase;
      tr.appendChild(tdPhase);
      const tdMol = document.createElement("td");
      tr.appendChild(tdMol);
      tbody.appendChild(tr);
      smiles ? renderMol2D(tdMol, smiles, inchikey) : (tdMol.textContent = "-");
    } else if (type === "bioactivity") {
      const [chemblId, mtype, evidence, inchikey, smiles] = row;
      const tdId = document.createElement("td");
      tdId.innerHTML =
        '<a href="' +
        chemblUrl(chemblId) +
        '" target="_blank" rel="noopener">' +
        esc(chemblId) +
        "</a>";
      tr.appendChild(tdId);
      const tdType = document.createElement("td");
      tdType.textContent = mtype || "-";
      tr.appendChild(tdType);
      const tdEv = document.createElement("td");
      tdEv.textContent = evidence || "-";
      tr.appendChild(tdEv);
      const tdMol = document.createElement("td");
      tr.appendChild(tdMol);
      tbody.appendChild(tr);
      smiles ? renderMol2D(tdMol, smiles, inchikey) : (tdMol.textContent = "-");
    }
  }
}

// =============================
// Populate Gene
// =============================
async function populateGene(gene) {
  if (!chemblData[gene]) return;
  currentGene = gene;
  const sections = [
    { key: "approved_drugs", id: "approvedSection", tbl: "approvedTable", title: "FDA-approved Drugs" },
    { key: "clinical_trials", id: "clinicalSection", tbl: "clinicalTable", title: "Clinical and Preclinical Trials" },
    { key: "bioactivity", id: "bioactivitySection", tbl: "bioactivityTable", title: "Bioactivity Evidence" },
  ];
  for (const { key, id, tbl, title } of sections) {
    const section = document.getElementById(id);
    const header = section.querySelector(".collapsible-header .title");
    const count = section.querySelector(".collapsible-header .count");
    const tableBody = document.getElementById(tbl);
    const data = chemblData[gene][key];
    if (!data || !data.length) {
      section.style.display = "none";
      continue;
    }
    header.textContent = title;
    if (key === "approved_drugs") {
      const grouped = groupApprovedRows(data);
      count.textContent = grouped.length;
      section.style.display = "block";
      await buildApprovedTable(tableBody, grouped);
    } else {
      count.textContent = data.length;
      section.style.display = "block";
      await buildTable(tableBody, data, key);
    }
  }
}

// =============================
// DOM Ready + Secure Theme Sync
// =============================
document.addEventListener("DOMContentLoaded", async () => {
  try {
    RDKit = await initRDKitModule();
  } catch (e) {
    console.warn("RDKit init failed:", e);
  }

  const geneSelect = document.getElementById("geneSelect");
  const validGenes = [];
  Object.keys(chemblData || {}).forEach((g) => {
    const d = chemblData[g];
    const hasData =
      (d?.approved_drugs?.length || 0) +
      (d?.clinical_trials?.length || 0) +
      (d?.bioactivity?.length || 0);
    if (hasData) {
      const opt = document.createElement("option");
      opt.value = g;
      opt.textContent = g;
      geneSelect.appendChild(opt);
      validGenes.push(g);
    }
  });

  if (validGenes.length > 0) {
    geneSelect.value = validGenes[0];
    populateGene(validGenes[0]);
  } else {
    geneSelect.style.display = "none";
    document.querySelectorAll(".collapsible-section").forEach((sec) => (sec.style.display = "none"));
  }

  document.querySelectorAll(".collapsible-header").forEach((header) => {
    const section = header.closest(".collapsible-section");
    const content = section.querySelector(".collapsible-content");
    const icon = header.querySelector(".toggle-icon");
    header.addEventListener("click", () => {
      const isOpen = content.classList.toggle("open");
      header.classList.toggle("active", isOpen);
      icon.textContent = isOpen ? "▾" : "▸";
      content.style.maxHeight = "";
    });
    content.classList.remove("open");
    content.style.maxHeight = "";
    icon.textContent = "▸";
  });

  geneSelect.addEventListener("change", (e) => populateGene(e.target.value));

  // Secure Theme Sync (no same-origin required)
  window.addEventListener("message", (e) => {
    if (!e.data || typeof e.data.theme !== "string") return;
    if (e.data.theme === "light") document.body.classList.add("light");
    else document.body.classList.remove("light");
  });
});
