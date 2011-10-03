WWW = /proj/web-cxc-dmz/htdocs/contrib/deproject

.PHONY: m87tar docs dist

dist:
	rm -rf dist
	python setup.py sdist

m87tar: 
	tar zcvf examples/m87.tar.gz examples/m87

docs:
	cd docs; \
	make html

install:
	rsync -av docs/.build/html/ $(WWW)/
	rsync -av examples/m87.tar.gz $(WWW)/downloads/
	cp dist/deproject-*.tar.gz $(WWW)/downloads/deproject.tar.gz
