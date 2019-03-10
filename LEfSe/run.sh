#!usr/bin/bash
format_input.py  lefse_input.tsv  lefse_input.in -c 1  -u 2 -o 1000000
run_lefse.py lefse_input.in  lefse.res
plot_res.py lefse.res  lefse.png --format png --dpi 300
plot_cladogram.py lefse.res lefse.cladogram.png --format png --dpi 300 --labeled_start_lev 2  --labeled_stop_lev 6 --abrv_start_lev 4 --abrv_stop_lev 6 
