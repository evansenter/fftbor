commands = Dir["*.txt"].map { |file| "time ~/Source/fftbor/FFTbor2D -R 75 -S -P 8 #{file} > #{File.basename(file, '.txt')}.out" }
run_pbs_job(commands, "giegerich_seqs_analysis_fftbor2d_R_75")
