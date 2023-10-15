#!/scinet/gpc/bin/make -Rf

# 14.12.2010 12:02:54 EST
# Harm van Bakel <hvbakel@gmail.com>

ifndef batch
$(error missing required value 'batch')
endif

all: $(addsuffix .lock,$(batch))

%.lock: ../0-mercounts/test.nmers.obt.fasta ../test.gkpStore
	./overlap.sh $(subst .lock,,$@)
