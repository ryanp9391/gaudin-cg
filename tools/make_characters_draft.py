import os
import shutil
import subprocess
import tempfile
import zipfile
import xml.etree.ElementTree as ET

SOURCE = "/home/paul/git_projects/Talks/Characters-JC.odp"
OUTPUT = "/home/paul/git_projects/Talks/Characters-JC-draft.odp"
BACKUP = "/home/paul/git_projects/Talks/Characters-JC-original.odp"
SVG_SOURCE = "/home/paul/git_projects/gaudin-CG/assets/coadjoint-orbit.svg"

NS = {
    "draw": "urn:oasis:names:tc:opendocument:xmlns:drawing:1.0",
    "office": "urn:oasis:names:tc:opendocument:xmlns:office:1.0",
    "text": "urn:oasis:names:tc:opendocument:xmlns:text:1.0",
    "svg": "urn:oasis:names:tc:opendocument:xmlns:svg-compatible:1.0",
    "xlink": "http://www.w3.org/1999/xlink",
}
for prefix, uri in NS.items():
    ET.register_namespace(prefix, uri)

slides = [
    ["CFT characters from localization", "Journal club", "Paul Ryan"],
    ["Based on arXiv:2609.05211", "Lorenz Eberhardt · Duarte Fragoso · Hessel Posthuma", "10–15 minute overview"],
    ["What is a character?", "Take a group or algebra G", "Choose a representation V", "Character: χV(g) = TrV(g)", "A state-counting function", "In CFT: a graded trace over states"],
    ["Why do characters matter?", "They package the spectrum of a representation", "They record how states transform under symmetries", "In 2d CFT, characters are partition functions", "The challenge: compute them efficiently"],
    ["The geometric idea", "Replace state counting by geometry", "coadjoint orbit  →  quantization  →  representation", "                         →  character", "The paper adapts the orbit method to CFT algebras"],
    ["Finite-dimensional warm-up: SU(2)", "A coadjoint orbit is a sphere", "S² ≅ CP¹ is a classical phase space", "Quantization gives H⁰(CP¹, O(n))", "This is the (n + 1)-dimensional representation", "Geometry already knows the representation theory"],
    ["Localization computes the trace", "Character = equivariant index", "χ = ∫orbit equivariant Todd class × line-bundle contribution", "Localization turns the integral into a sum over fixed points", "Only local fixed-point data is needed"],
    ["A familiar answer reappears", "For compact groups, localization reproduces Weyl’s character formula", "Fixed points ↔ Weyl-group data", "Local weights ↔ the denominator", "This is the finite-dimensional proof of concept"],
    ["From groups to CFT algebras", "Virasoro and affine Kac–Moody orbits are infinite-dimensional", "The Kähler geometry is formal rather than fully rigorous", "But the same logical chain survives", "orbit  →  quantization  →  localization  →  character"],
    ["Example: Virasoro characters", "The orbit method reproduces the character of a Virasoro module", "Fixed-point contributions generate the oscillator product", "Null-state information is encoded geometrically", "The result matches the standard CFT character"],
    ["How far does it extend?", "Affine Kac–Moody algebras", "N = 1, 2 and small N = 4 super-Virasoro algebras", "Fermionic directions lead to a super-index", "The same method gives a uniform derivation"],
    ["Take-home message", "Characters are traces that encode representation spectra", "Coadjoint orbits provide geometric models of representations", "Equivariant localization reduces character computations to fixed points", "The method unifies finite-dimensional representation theory and CFT", "Questions?"],
]


def add_text_shape(page, style, text_style, x, y, width, height, lines, name):
    shape = ET.Element("{%s}custom-shape" % NS["draw"], {
        "{%s}name" % NS["draw"]: name,
        "{%s}style-name" % NS["draw"]: style,
        "{%s}text-style-name" % NS["draw"]: text_style,
        "{%s}layer" % NS["draw"]: "layout",
        "{%s}width" % NS["svg"]: width,
        "{%s}height" % NS["svg"]: height,
        "{%s}x" % NS["svg"]: x,
        "{%s}y" % NS["svg"]: y,
    })
    geometry = ET.SubElement(shape, "{%s}enhanced-geometry" % NS["draw"], {
        "{%s}mirror-horizontal" % NS["draw"]: "false",
        "{%s}mirror-vertical" % NS["draw"]: "false",
        "{%s}viewBox" % NS["svg"]: "0 0 0 0",
        "{%s}text-areas" % NS["draw"]: "0 0 ?f3 ?f2",
        "{%s}type" % NS["draw"]: "ooxml-rect",
        "{%s}enhanced-path" % NS["draw"]: "M 0 0 L ?f3 0 ?f3 ?f2 0 ?f2 Z N",
    })
    ET.SubElement(geometry, "{%s}equation" % NS["draw"], {"{%s}name" % NS["draw"]: "f2", "{%s}formula" % NS["draw"]: "logheight"})
    ET.SubElement(geometry, "{%s}equation" % NS["draw"], {"{%s}name" % NS["draw"]: "f3", "{%s}formula" % NS["draw"]: "logwidth"})
    for line in lines:
        paragraph = ET.SubElement(shape, "{%s}p" % NS["text"], {"{%s}style-name" % NS["text"]: "P3" if style == "gr12" else "P10"})
        span = ET.SubElement(paragraph, "{%s}span" % NS["text"], {"{%s}style-name" % NS["text"]: "T2" if style == "gr12" else "T6"})
        span.text = line
    page.append(shape)


def add_image(page):
    frame = ET.Element("{%s}frame" % NS["draw"], {
        "{%s}style-name" % NS["draw"]: "gr16",
        "{%s}name" % NS["draw"]: "Coadjoint orbit diagram",
        "{%s}layer" % NS["draw"]: "layout",
        "{%s}width" % NS["svg"]: "12.5cm",
        "{%s}height" % NS["svg"]: "4.5cm",
        "{%s}x" % NS["svg"]: "12.0cm",
        "{%s}y" % NS["svg"]: "3.3cm",
    })
    image = ET.SubElement(frame, "{%s}image" % NS["draw"], {
        "{%s}href" % NS["xlink"]: "Pictures/coadjoint-orbit.svg",
        "{%s}type" % NS["xlink"]: "simple",
        "{%s}show" % NS["xlink"]: "embed",
        "{%s}actuate" % NS["xlink"]: "onLoad",
    })
    page.append(frame)

    work = tempfile.mkdtemp(prefix="characters-odp-")
    with zipfile.ZipFile(SOURCE) as archive:
        archive.extractall(work)
    shutil.copy2(SVG_SOURCE, os.path.join(work, "Pictures", "coadjoint-orbit.svg"))

    content_path = os.path.join(work, "content.xml")
    tree = ET.parse(content_path)
    root = tree.getroot()
    presentation = root.find(".//office:presentation", NS)
    pages = presentation.findall("draw:page", NS)

    # Keep the original title slide's composition, but use clean, uniform
    # text boxes on the content slides so diagram-specific source objects do
    # not interfere with the new narrative.
    title_page = pages[0]
    paragraphs = title_page.findall(".//text:p", NS)
    nonempty = [p for p in paragraphs if "".join(p.itertext()).strip()]
    for p in paragraphs:
        for child in list(p):
            p.remove(child)
        p.text = ""
    for p, value in zip(nonempty, slides[0]):
        p.text = value

    for index, (page, replacement) in enumerate(zip(pages[1: len(slides)], slides[1:]), 2):
        for child in list(page):
            if child.tag != "{%s}forms" % "urn:oasis:names:tc:opendocument:xmlns:form:1.0":
                page.remove(child)
        add_text_shape(page, "gr12", "P4", "1.5cm", "0.8cm", "23.5cm", "1.4cm", [replacement[0]], "Draft title %d" % index)
        body = replacement[1:]
        if index == 3:
            add_text_shape(page, "gr15", "P11", "3.2cm", "3.0cm", "20.8cm", "7.0cm", body, "Definition box")
        elif index == 5:
            add_text_shape(page, "gr16", "P11", "2.0cm", "2.8cm", "9.0cm", "8.0cm", body, "Draft body %d" % index)
            add_image(page)
        else:
            add_text_shape(page, "gr16", "P11", "2.0cm", "2.8cm", "22.5cm", "9.0cm", body, "Draft body %d" % index)

    manifest_path = os.path.join(work, "META-INF", "manifest.xml")
    manifest_tree = ET.parse(manifest_path)
    manifest_root = manifest_tree.getroot()
    manifest_ns = "urn:oasis:names:tc:opendocument:xmlns:manifest:1.0"
    ET.register_namespace("manifest", manifest_ns)
    ET.SubElement(manifest_root, "{%s}file-entry" % manifest_ns, {
        "{%s}full-path" % manifest_ns: "Pictures/coadjoint-orbit.svg",
        "{%s}media-type" % manifest_ns: "image/svg+xml",
    })
    manifest_tree.write(manifest_path, encoding="UTF-8", xml_declaration=True)

    for page in pages[len(slides) :]:
        presentation.remove(page)

    tree.write(content_path, encoding="UTF-8", xml_declaration=True)
    if not os.path.exists(BACKUP):
        shutil.copy2(SOURCE, BACKUP)

    temp_output = OUTPUT + ".tmp"
    if os.path.exists(temp_output):
        os.remove(temp_output)
    with zipfile.ZipFile(temp_output, "w") as archive:
        for dirname, _, filenames in os.walk(work):
            for filename in filenames:
                path = os.path.join(dirname, filename)
                arcname = os.path.relpath(path, work)
                compression = zipfile.ZIP_STORED if arcname == "mimetype" else zipfile.ZIP_DEFLATED
                archive.write(path, arcname, compress_type=compression)
    os.replace(temp_output, OUTPUT)
    shutil.rmtree(work, ignore_errors=True)

print(OUTPUT)
