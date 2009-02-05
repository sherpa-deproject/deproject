WWW = /proj/web-cxc/htdocs/contrib/deproject

.PHONY: help clean html web pickle htmlhelp latex changes linkcheck docs

install:
	rsync -av docs/.build/html/ $(WWW)/
	rsync -av dist/deproject-*.tar.gz examples/m87.tar.gz $(WWW)/downloads/

m87tar: 
	tar zcvf m87.tar.gz examples/m87

docs:
	cd docs; \
	make html
