.PHONY: purge clean check

default:
	latexmk -f -pdf

purge: clean
	latexmk -C
	$(RM) -r out

clean:
	latexmk -c -bibtex
	$(RM) *.synctex.gz
