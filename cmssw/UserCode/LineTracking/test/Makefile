tags = simu13TeV

all :

sub :
	@ for tag in $(tags) ; do \
	 rm -rf $$tag ; \
	 crab -cfg $$tag.cfg -create ; \
	 crab -submit -c $$tag ; \
	done

stat :
	@ for tag in $(tags) ; do \
	 crab -status -c $$tag ; \
	done

get :
	@ for tag in $(tags) ; do \
	 crab -get -c $$tag ; \
	done

kill :
	@ for tag in $(tags) ; do \
	 crab -kill all -c $$tag ; \
	done
