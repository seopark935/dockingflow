// Convert the DockingFlow manual (a Markdown subset) into a formatted .docx.
//   node md2docx.js <input.md> <output.docx>
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, BorderStyle, AlignmentType, TableOfContents, LevelFormat,
  PageBreak, ExternalHyperlink,
} = require("docx");

const [, , IN, OUT] = process.argv;
const SRC = fs.readFileSync(IN, "utf8").split("\n");

const BODY = "Calibri", MONO = "Consolas";
const ACCENT = "1F4E79", MUTED = "595959", CODEBG = "F2F3F5", HEADBG = "DCE6F1";
const PAGE_W = 12240, MARGIN = 1080;          // US Letter, 0.75" margins
const CONTENT_W = PAGE_W - 2 * MARGIN;

// ---------- inline markdown: **bold**, `code`, [text](url), *italic* ----------
function inline(text, base = {}) {
  const runs = [];
  const re = /(\*\*[^*]+\*\*|`[^`]+`|\[[^\]]+\]\([^)]+\)|\*[^*]+\*)/g;
  let last = 0, m;
  const push = (t, extra) => { if (t) runs.push(new TextRun({ text: t, font: BODY, ...base, ...extra })); };
  while ((m = re.exec(text)) !== null) {
    push(clean(text.slice(last, m.index)));
    const tok = m[0];
    if (tok.startsWith("**")) push(clean(tok.slice(2, -2)), { bold: true });
    else if (tok.startsWith("`")) runs.push(new TextRun({
      text: tok.slice(1, -1), font: MONO, size: 19, shading: { type: ShadingType.CLEAR, fill: CODEBG }, ...base }));
    else if (tok.startsWith("[")) {
      const mm = /\[([^\]]+)\]\(([^)]+)\)/.exec(tok);
      if (/^https?:/.test(mm[2])) {
        runs.push(new ExternalHyperlink({
          link: mm[2],
          children: [new TextRun({ text: clean(mm[1]), font: BODY, color: "0563C1", underline: {}, ...base })],
        }));
      } else push(clean(mm[1]));           // internal anchors: keep the words, drop the link
    } else push(clean(tok.slice(1, -1)), { italics: true });
    last = m.index + tok.length;
  }
  push(clean(text.slice(last)));
  return runs.length ? runs : [new TextRun({ text: "", font: BODY })];
}

// Markdown escapes and the few non-ASCII glyphs Word renders oddly in Consolas.
const clean = (t) => t.replace(/\\([|*_`])/g, "$1");

const para = (text, opts = {}) => new Paragraph({
  children: inline(text), spacing: { after: 120, line: 276 }, ...opts,
});

function codeBlock(lines) {
  // One paragraph per line, shaded, kept together.
  return lines.map((ln, i) => new Paragraph({
    children: [new TextRun({ text: ln || " ", font: MONO, size: 17 })],
    shading: { type: ShadingType.CLEAR, fill: CODEBG },
    spacing: { after: i === lines.length - 1 ? 160 : 0, before: i === 0 ? 60 : 0, line: 240 },
    keepNext: i < lines.length - 1,
    indent: { left: 180, right: 180 },
  }));
}

function table(rows) {
  const cols = rows[0].length;
  const widths = Array(cols).fill(Math.floor(CONTENT_W / cols));
  widths[cols - 1] = CONTENT_W - widths[0] * (cols - 1);   // make them sum exactly
  return new Table({
    columnWidths: widths,
    width: { size: CONTENT_W, type: WidthType.DXA },
    rows: rows.map((cells, r) => new TableRow({
      tableHeader: r === 0,
      children: cells.map((cell, c) => new TableCell({
        width: { size: widths[c], type: WidthType.DXA },
        shading: r === 0 ? { type: ShadingType.CLEAR, fill: HEADBG } : undefined,
        margins: { top: 60, bottom: 60, left: 100, right: 100 },
        children: [new Paragraph({
          children: inline(cell, r === 0 ? { bold: true } : {}),
          spacing: { after: 0, line: 252 },
        })],
      })),
    })),
  });
}

// ---------- block parsing ----------
const out = [];
let i = 0;
const bullets = [], numbers = [];

while (i < SRC.length) {
  const line = SRC[i];

  if (/^\s*$/.test(line)) { i++; continue; }

  // horizontal rule -> page break (the manual uses --- between major sections)
  if (/^---+\s*$/.test(line)) {
    out.push(new Paragraph({ children: [new PageBreak()] }));
    i++; continue;
  }

  // fenced code
  if (/^```/.test(line)) {
    const body = [];
    i++;
    while (i < SRC.length && !/^```/.test(SRC[i])) body.push(SRC[i++]);
    i++;
    out.push(...codeBlock(body));
    continue;
  }

  // table
  if (/^\|/.test(line) && /^\|[\s:|-]+\|?\s*$/.test(SRC[i + 1] || "")) {
    const rows = [];
    const cells = (ln) => ln.replace(/^\||\|$/g, "").split(/(?<!\\)\|/).map((c) => c.trim());
    rows.push(cells(line));
    i += 2;
    while (i < SRC.length && /^\|/.test(SRC[i])) rows.push(cells(SRC[i++]));
    out.push(table(rows));
    out.push(new Paragraph({ text: "", spacing: { after: 120 } }));
    continue;
  }

  // headings
  const h = /^(#{1,6})\s+(.*)$/.exec(line);
  if (h) {
    const level = h[1].length;
    // Headings keep their own numbering ("2.1 Connect ...") verbatim: stripping
    // it mangles multi-level numbers.
    out.push(new Paragraph({
      children: inline(h[2]),
      heading: level === 1 ? HeadingLevel.TITLE
        : level === 2 ? HeadingLevel.HEADING_1
        : level === 3 ? HeadingLevel.HEADING_2 : HeadingLevel.HEADING_3,
      spacing: { before: level <= 2 ? 300 : 240, after: 120 },
      keepNext: true,
    }));
    i++; continue;
  }

  // blockquote
  if (/^>\s?/.test(line)) {
    const body = [];
    while (i < SRC.length && /^>\s?/.test(SRC[i])) body.push(SRC[i++].replace(/^>\s?/, ""));
    out.push(new Paragraph({
      children: inline(body.join(" ")),
      indent: { left: 360 },
      border: { left: { style: BorderStyle.SINGLE, size: 12, color: ACCENT, space: 12 } },
      spacing: { after: 160, before: 80, line: 276 },
    }));
    continue;
  }

  // screenshot placeholder
  const img = /^!\[([^\]]*)\]\(([^)]+)\)/.exec(line);
  if (img) {
    out.push(new Paragraph({
      children: [new TextRun({ text: `[ Screenshot: ${img[1]} — insert ${img[2]} ]`, font: BODY, italics: true, color: MUTED })],
      alignment: AlignmentType.CENTER,
      shading: { type: ShadingType.CLEAR, fill: CODEBG },
      spacing: { before: 120, after: 160 },
      border: {
        top: { style: BorderStyle.DASHED, size: 6, color: "AAAAAA", space: 8 },
        bottom: { style: BorderStyle.DASHED, size: 6, color: "AAAAAA", space: 8 },
        left: { style: BorderStyle.DASHED, size: 6, color: "AAAAAA", space: 8 },
        right: { style: BorderStyle.DASHED, size: 6, color: "AAAAAA", space: 8 },
      },
    }));
    i++; continue;
  }

  // lists
  const b = /^(\s*)[-*]\s+(.*)$/.exec(line);
  const n = /^(\s*)(\d+)\.\s+(.*)$/.exec(line);
  if (b || n) {
    const indentLevel = Math.min(Math.floor(((b || n)[1] || "").length / 2), 2);
    let text = b ? b[2] : n[3];
    // fold continuation lines into the same item
    while (i + 1 < SRC.length && /^\s{2,}\S/.test(SRC[i + 1]) && !/^\s*([-*]|\d+\.)\s/.test(SRC[i + 1])
           && !/^\s*```/.test(SRC[i + 1])) {
      text += " " + SRC[++i].trim();
    }
    out.push(new Paragraph({
      children: inline(text),
      numbering: { reference: b ? "bullets" : "numbers", level: indentLevel },
      spacing: { after: 80, line: 276 },
    }));
    i++; continue;
  }

  // plain paragraph (join wrapped lines)
  let text = line.trim();
  while (i + 1 < SRC.length && SRC[i + 1].trim() && !/^(#{1,6}\s|```|\||>|!\[|\s*[-*]\s|\s*\d+\.\s|---+\s*$)/.test(SRC[i + 1])) {
    text += " " + SRC[++i].trim();
  }
  out.push(para(text));
  i++;
}

// ---------- document ----------
const headingStyle = (size) => ({
  run: { font: BODY, size, bold: true, color: ACCENT },
  paragraph: { spacing: { before: 280, after: 120 }, keepNext: true },
});

const doc = new Document({
  creator: "DockingFlow",
  title: "DockingFlow Lab Manual",
  description: "Running structure-based virtual screens with DockingFlow",
  styles: {
    default: {
      document: { run: { font: BODY, size: 21 }, paragraph: { spacing: { after: 120, line: 276 } } },
      title: {
        run: { font: BODY, size: 52, bold: true, color: ACCENT },
        paragraph: { spacing: { before: 0, after: 240 }, keepNext: true },
      },
      heading1: headingStyle(32),
      heading2: headingStyle(26),
      heading3: headingStyle(23),
    },
  },
  numbering: {
    config: [
      {
        reference: "bullets",
        levels: [0, 1, 2].map((lvl) => ({
          level: lvl, format: LevelFormat.BULLET, text: ["•", "◦", "▪"][lvl],
          alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 360 + lvl * 360, hanging: 260 } } },
        })),
      },
      {
        reference: "numbers",
        levels: [0, 1, 2].map((lvl) => ({
          level: lvl, format: LevelFormat.DECIMAL, text: `%${lvl + 1}.`,
          alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 360 + lvl * 360, hanging: 260 } } },
        })),
      },
    ],
  },
  sections: [{
    properties: { page: { size: { width: PAGE_W, height: 15840 }, margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN } } },
    children: [
      new Paragraph({
        children: [new TextRun({ text: "Contents", font: BODY, size: 32, bold: true, color: ACCENT })],
        spacing: { before: 0, after: 160 }, keepNext: true,
      }),
      new TableOfContents("Contents", { hyperlink: true, headingStyleRange: "1-3" }),
      new Paragraph({ children: [new PageBreak()] }),
      ...out,
    ],
  }],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync(OUT, buf);
  console.log(`wrote ${OUT} (${(buf.length / 1024).toFixed(0)} KB, ${out.length} blocks)`);
});
