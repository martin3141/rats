main: generate_results.R gsh_mpress_data.Rda
	Rscript	generate_results.R

clean:
	rm *.eps
