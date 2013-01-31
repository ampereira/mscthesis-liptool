CMD=latexmk
FLAGS=-f -pdf -quiet

default:
	$(CMD) $(FLAGS)

once:


%: tex/%.tex
	$(CMD) $(FLAGS) tex/$@.tex

clean:
	latexmk -C
