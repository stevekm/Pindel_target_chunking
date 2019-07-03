SHELL:=/bin/bash
export NXF_VER:=19.01.0
./nextflow:
	curl -fsSL get.nextflow.io | bash
install: ./nextflow

# ~~~~~ RUN PIPELINE ~~~~~ #
TIMESTAMP:=$(shell date +%s)
TIMESTAMP_str:=$(shell date +"%Y-%m-%d-%H-%M-%S")
ABSDIR:=$(shell python -c 'import os; print(os.path.realpath("."))')
DIRNAME:=$(shell python -c 'import os; print(os.path.basename(os.path.realpath(".")))')
LOGDIR:=logs
LOGDIRABS:=$(shell python -c 'import os; print(os.path.realpath("$(LOGDIR)"))')
LOGID:=$(TIMESTAMP)
LOGFILEBASE:=log.$(LOGID).out
LOGFILE:=$(LOGDIR)/$(LOGFILEBASE)
# start run with logging
run: install
	@log_file="$(LOGDIR)/nextflow.$(LOGID).stdout.log" ; \
	echo ">>> Running with stdout log file: $(LOGFILE)" ; \
	$(MAKE) run-recurse 2>&1 | tee -a "$(LOGFILE)" ; \
	echo ">>> Run completed at $$(date +"%Y-%m-%d %H:%M:%S"), stdout log file: $(LOGFILE)"

# try to automatically determine which 'run' recipe to use based on hostname
run-recurse:
	./nextflow run main.nf -resume $(EP)

# submit the parent Nextflow process to Big Purple HPC as a cluster job
SUBJOBNAME:=$(DIRNAME)
SUBLOG:=$(LOGDIRABS)/slurm-%j.$(LOGFILEBASE)
SUBQ:=intellispace
SUBTIME:=--time=5-00:00:00
SUBTHREADS:=4
SUBEP:=
NXF_NODEFILE:=.nextflow.node
NXF_JOBFILE:=.nextflow.jobid
NXF_PIDFILE:=.nextflow.pid
NXF_SUBMIT:=.nextflow.submitted
NXF_SUBMITLOG:=.nextflow.submitted.log
REMOTE:=
PID:=
submit: install
	@if [ -e "$(NXF_SUBMIT)" ]; then echo ">>> ERROR: An instance of the pipeline has already been submitted"; exit 1 ; \
	else echo ">>> Running submit-bigpurple"; $(MAKE) submit-bigpurple ; fi

submit-bigpurple:
	@touch "$(NXF_SUBMIT)" && \
	printf "#!/bin/bash\n \
	make submit-bigpurple-run TIMESTAMP=$(TIMESTAMP) $(SUBEP)" | \
	sbatch -D "$(ABSDIR)" -o "$(SUBLOG)" -J "$(SUBJOBNAME)" -p "$(SUBQ)" $(SUBTIME) --ntasks-per-node=1 -c "$(SUBTHREADS)" --export=HOSTNAME /dev/stdin | tee >(sed 's|[^[:digit:]]*\([[:digit:]]*\).*|\1|' > '$(NXF_JOBFILE)')

submit-bigpurple-run:
	if [ -e "$(NXF_NODEFILE)" -a -e "$(NXF_PIDFILE)" ]; then paste "$(NXF_NODEFILE)" "$(NXF_PIDFILE)" >> $(NXF_SUBMITLOG); fi ; \
	echo "$${SLURMD_NODENAME}" > "$(NXF_NODEFILE)" && \
	$(MAKE) run HOSTNAME="bigpurple" LOGID="$(TIMESTAMP)" EP='-bg' && \
	if [ -e "$(NXF_SUBMIT)" ]; then rm -f "$(NXF_SUBMIT)"; fi

# kill the Nextflow process running inside the submitted HPC job, in order to cleanly end the pipeline
kill: PID=$(shell head -1 "$(NXF_PIDFILE)")
kill: REMOTE=$(shell head -1 "$(NXF_NODEFILE)")
kill: $(NXF_NODEFILE) $(NXF_PIDFILE)
	ssh "$(REMOTE)" 'kill $(PID)'


PRE:=
RECDIR:=recorded-runs/$(PRE)$(DIRNAME)_$(TIMESTAMP_str)
STDOUTLOGPATH:=
STDOUTLOG:=
ALL_LOGS:=
# save a hard copy of the latest pipeline execution logs
record: STDOUTLOGPATH=$(shell ls -d -1t $(LOGDIR)/log.*.out | head -1 | python -c 'import sys, os; print(os.path.realpath(sys.stdin.readlines()[0].strip()))' )
record: STDOUTLOG=$(shell basename "$(STDOUTLOGPATH)")
record: ALL_LOGS=$(shell find "$(LOGDIR)" -type f -name '*$(STDOUTLOG)*')
record:
	@mkdir -p "$(RECDIR)" && \
	cp -a *.html trace.txt .nextflow.log main.nf nextflow.config "$(RECDIR)/" && \
	for item in $(ALL_LOGS); do cp -a "$${item}" "$(RECDIR)/"; done ; \
	echo ">>> Copied execution reports and logs to: $(RECDIR)"

# ~~~~~ CLEANUP ~~~~~ #
# commands to clean out items in the current directory after running the pipeline
clean-traces:
	rm -f trace*.txt.*

clean-logs:
	rm -f .nextflow.log.*

clean-reports:
	rm -f *.html.*

clean-flowcharts:
	rm -f *.dot.*

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

clean-old-stdout-logs:
	find . -maxdepth 1 -mindepth 1 -type f -name "nextflow.*.stdout.log" | sort -r | tail -n +2 | xargs rm -f

clean-old-failed-logs:
	find . -maxdepth 1 -mindepth 1 -type f -name "failed.*.tsv" | sort -r | tail -n +2 | xargs rm -f

# clean all files produced by previous pipeline runs
clean: clean-logs clean-traces clean-reports clean-flowcharts clean-old-stdout-logs clean-old-failed-logs

# clean all files produced by all pipeline runs
clean-all: clean clean-output clean-work
	[ -d .nextflow ] && mv .nextflow .nextflowold && rm -rf .nextflowold &
	rm -f .nextflow.log
	rm -f *.png
	rm -f trace*.txt*
	rm -f *.html*
	rm -f flowchart*.dot
	rm -f nextflow.*.stdout.log
