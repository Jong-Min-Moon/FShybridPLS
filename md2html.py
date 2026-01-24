
import markdown
import sys
import os

if len(sys.argv) < 3:
    print("Usage: python3 md2html.py <input.md> <output.html>")
    sys.exit(1)

input_path = sys.argv[1]
output_path = sys.argv[2]

with open(input_path, 'r', encoding='utf-8') as f:
    text = f.read()

# Convert markdown to HTML
# 'extra' enables tables, definition lists, etc.
# 'codehilite' enables code highlighting
try:
    html_content = markdown.markdown(text, extensions=['extra', 'codehilite', 'toc'])
except ImportError:
    # Fallback if extensions are missing, though 'extra' usually comes with standard install
    html_content = markdown.markdown(text)

# rudimentary mathjax or simple replacements could be added here if needed, 
# but for PDF via cupsfilter, simple HTML is best.
# Latex math might not render in cupsfilter without JS.
# We will assume the user accepts the latex source code representation or we attempt a simple replacement.
# Since I fixed the LaTeX to be standard, it will show as $E=mc^2$ in the PDF which is acceptable for code documentation often.

html_full = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
body {{ 
    font-family: Helvetica, sans-serif; 
    font-size: 11pt; 
    line-height: 1.6; 
    max-width: 800px;
    margin: 0 auto;
    padding: 20px;
}}
pre {{ 
    background: #f4f4f4; 
    padding: 10px; 
    border: 1px solid #ddd;
    border-radius: 5px; 
    white-space: pre-wrap; 
    overflow-x: auto;
}}
code {{ 
    font-family: monospace; 
    background: #f4f4f4; 
    padding: 2px 4px; 
    border-radius: 3px; 
    font-size: 0.9em;
}}
h1, h2, h3, h4, h5, h6 {{ 
    color: #333; 
    margin-top: 1.5em;
    margin-bottom: 0.5em;
}}
h1 {{ border-bottom: 2px solid #eee; padding-bottom: 10px; }}
h2 {{ border-bottom: 1px solid #eee; padding-bottom: 5px; }}
table {{ 
    border-collapse: collapse; 
    width: 100%; 
    margin-bottom: 1em; 
}}
th, td {{ 
    border: 1px solid #ddd; 
    padding: 8px; 
    text-align: left; 
}}
th {{ background-color: #f7f7f7; font-weight: bold; }}
blockquote {{
    border-left: 4px solid #ddd;
    padding-left: 10px;
    color: #666;
    margin-left: 0;
}}
img {{ max-width: 100%; }}
</style>
</head>
<body>
{html_content}
</body>
</html>
"""

with open(output_path, 'w', encoding='utf-8') as f:
    f.write(html_full)
