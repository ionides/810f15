%.html: %.Rmd
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

