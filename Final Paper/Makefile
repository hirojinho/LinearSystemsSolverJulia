# Makefile for LaTeX Project - DMD Final Paper
# MAT-55 - Álgebra Linear Computacional

# Variables
MAIN = DMD_Trabalho_Final_MAT55
TEX_FILES = $(wildcard *.tex)
PDF_FILES = $(MAIN).pdf
LATEX_CMD = pdflatex
BIBTEX_CMD = bibtex

# Default target
all: $(PDF_FILES)

# Main compilation rule
$(MAIN).pdf: $(MAIN).tex $(TEX_FILES)
	@echo "Compiling LaTeX document..."
	$(LATEX_CMD) $(MAIN).tex
	@echo "Running second pass for cross-references..."
	$(LATEX_CMD) $(MAIN).tex
	@echo "Compilation complete!"

# Alternative compilation with bibliography (if needed)
full: $(MAIN).tex $(TEX_FILES)
	@echo "Full compilation with bibliography..."
	$(LATEX_CMD) $(MAIN).tex
	$(BIBTEX_CMD) $(MAIN)
	$(LATEX_CMD) $(MAIN).tex
	$(LATEX_CMD) $(MAIN).tex
	@echo "Full compilation complete!"

# Quick compilation (single pass)
quick: $(MAIN).tex
	@echo "Quick compilation (single pass)..."
	$(LATEX_CMD) $(MAIN).tex

# Clean temporary files
clean:
	@echo "Cleaning temporary files..."
	rm -f *.aux *.log *.out *.toc *.lot *.lof *.bbl *.blg *.fls *.fdb_latexmk *.synctex.gz *.nav *.snm *.vrb

# Clean all generated files including PDF
distclean: clean
	@echo "Cleaning all generated files..."
	rm -f $(PDF_FILES)

# View PDF
view: $(PDF_FILES)
	@echo "Opening PDF..."
	open $(MAIN).pdf

# Continuous compilation (watches for changes)
watch:
	@echo "Watching for changes... (Press Ctrl+C to stop)"
	while true; do \
		inotifywait -e modify $(TEX_FILES) 2>/dev/null || sleep 2; \
		make quick; \
	done

# Check LaTeX installation
check:
	@echo "Checking LaTeX installation..."
	@which pdflatex >/dev/null && echo "✓ pdflatex found" || echo "✗ pdflatex not found"
	@which bibtex >/dev/null && echo "✓ bibtex found" || echo "✗ bibtex not found"

# Word count (approximate)
wordcount: $(MAIN).tex
	@echo "Approximate word count:"
	@detex $(MAIN).tex | wc -w

# Spell check (if aspell is available)
spell: $(MAIN).tex
	@echo "Checking spelling..."
	@which aspell >/dev/null && aspell check --lang=pt_BR $(MAIN).tex || echo "aspell not available"

# Archive project
archive:
	@echo "Creating project archive..."
	@DATE=$$(date +%Y%m%d_%H%M%S); \
	tar -czf ../DMD_Paper_$$DATE.tar.gz --exclude-from=.gitignore .
	@echo "Archive created in parent directory"

# Help
help:
	@echo "Available targets:"
	@echo "  all       - Compile PDF (default)"
	@echo "  full      - Full compilation with bibliography"
	@echo "  quick     - Quick compilation (single pass)"
	@echo "  clean     - Remove temporary files"
	@echo "  distclean - Remove all generated files"
	@echo "  view      - Open PDF viewer"
	@echo "  watch     - Continuous compilation"
	@echo "  check     - Check LaTeX installation"
	@echo "  wordcount - Count words in document"
	@echo "  spell     - Spell check (requires aspell)"
	@echo "  archive   - Create project archive"
	@echo "  help      - Show this help"

# Declare phony targets
.PHONY: all full quick clean distclean view watch check wordcount spell archive help

# Set default shell
SHELL := /bin/bash 