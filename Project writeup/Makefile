# Build the paper PDF from main.tex.  Two pdflatex passes to resolve refs.
# (No bibtex needed: bibliography is inline via \begin{thebibliography}.)

TEX = pdflatex -interaction=nonstopmode -halt-on-error

main.pdf: main.tex
	$(TEX) main.tex
	$(TEX) main.tex

watch:
	@echo "Re-builds on every save (requires inotifywait)."
	@while inotifywait -e modify main.tex 2>/dev/null; do \
	  $(MAKE) -s main.pdf; \
	done

clean:
	rm -f main.aux main.log main.out main.toc main.bbl main.blg main.fdb_latexmk main.fls main.synctex.gz

distclean: clean
	rm -f main.pdf

.PHONY: watch clean distclean
