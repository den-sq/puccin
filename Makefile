.PHONY: all
all: clean vendor style test


.PHONY: help
help:
	@echo 'Management commands for phylo:'
	@echo
	@echo 'Usage:'
	@echo '    make clean           Clean the directory tree.'
	@echo '    make test            Run tests on the project.'
	@echo '    make test/benchmark  Run benchmark tests on the project.'
	@echo '    make vendor          ensures dependencies are installed.'
	@echo

##############################################################################
# The following targets are used for aiding in development and CI for the 
# phylo source code
##############################################################################
.PHONY: clean
clean:
	cargo clean

.PHONY: style
style:
	rustfmt --version
    cargo fmt -- --write-mode=diff

.PHONY: test
test:
	cargo test

.PHONY: test/benchmark
test/benchmark:
	cargo bench

.PHONY: vendor
vendor:
	cargo update

