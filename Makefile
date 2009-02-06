WWW = /proj/web-cxc/htdocs/contrib/deproject

.PHONY: m87tar docs dist

dist:
	python setup.py sdist

m87tar: 
	tar zcvf examples/m87.tar.gz examples/m87

docs:
	cd docs; \
	make html

install:
	rsync -av docs/.build/html/ $(WWW)/
	rsync -av dist/deproject-*.tar.gz examples/m87.tar.gz $(WWW)/downloads/
