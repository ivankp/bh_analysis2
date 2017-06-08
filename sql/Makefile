DB := test.db
DZ := $(DB).gz

ADD := $(wildcard add_ntuples_*.py)

.PHONY: all clean
.PRECIOUS: $(DB)

all: $(DB)

$(DB):
	@sqlite3 $@ < init.sql
	@sqlite3 $@ <<< '.schema'
	@for add in $(ADD); do python $$add $@; done

%.gz: %
	@gzip $<

clean:
	@rm -vf $(DB)
