all: doc install check site

site:
	Rscript -e 'pkgdown::build_site()'

doc:
	Rscript -e 'devtools::document()'

check:
	Rscript -e 'devtools::check()'

build:
	Rscript -e 'devtools::build()'

install:
	Rscript -e 'devtools::install()'

README.md: README.Rmd
	Rscript -e 'rmarkdown::render("README.Rmd", output_format = "md_document")'

