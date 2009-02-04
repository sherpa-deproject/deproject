WWW = /proj/web-cxc/htdocs/contrib/deproject

.PHONY: help clean html web pickle htmlhelp latex changes linkcheck docs

wwwdocs: docs
	rsync -av docs/.build/html/ $(WWW)/

docs:
	cd docs; \
	make html
