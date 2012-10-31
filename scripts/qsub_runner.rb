#! /usr/bin/ruby
# Exampe commandfile:
# ./runRNAmutantsModified.py BC021677.1_1965-1834 GUGGAGGUCUAUUCCAAUGGGGCUUUUCCUGUAGCUGCAUGUUGUUGGAAACUCCUCAUAGACUAACUCUGUGGUUUUGCUUUACUCACAGGACUAUUAGUUAGGUCUGUGGGAAGGAACUACAAGACAGUU /cluster/home/dingyc/2008summer/tmp/RNAmutants/RNAmutantsOutput/RF00425/BC021677.1_1965-1834/BC021677.1_1965-1834

require "vienna_rna"

def run_job(action, prefix)
  content = <<-SH
      #!/bin/sh
      #PBS -l nodes=1:clotelab
      #PBS -o #{prefix}.log
      #PBS -e #{prefix}.err
      #PBS -q normal
      #PBS -l walltime=2400:00:00
      cd $PBS_O_WORKDIR
      #{action}
    SH
  
  File.open("PBS_#{prefix}.sh", "w") do |file|
    file.write(content.gsub(/^\s*/, ""))
  end
      
  puts "Submitting job: #{action}"
      
  %x|qsub "PBS_#{prefix}.sh"|
end

unless ARGV.length == 2
  puts "Usage: ./qsub_runner.rb [output_prefix] [command_file]"
else
  puts "Add a command to (optionally) clean everything up first."
  
  action = File.read(ARGV[-1]).split(/\n/).first
  prefix = ARGV[-2].gsub(/\s+/, "_")
  
  if !File.exist?("#{prefix}.log") && !File.exist?("#{prefix}.err")
    run_job(action, prefix)
  elsif File.exist?("#{prefix}.log") || File.exist?("#{prefix}.err")
    puts "Warning: #{prefix}.log or #{prefix}.err already exists. Continue? Y/N"
    
    if STDIN.gets =~ /y/i
      run_job(action, prefix)
    end
  end
end

# --------------------------------------------+
# RNAbor / FFTbor time benchmarking functions |
# --------------------------------------------+
def make_command_file(algorithm, size, &block)
  content = <<-SH
    #!/bin/sh
    #PBS -l nodes=1:clotelab
    #PBS -o #{size}_#{algorithm}.log
    #PBS -e #{size}_#{algorithm}.err
    #PBS -q normal
    #PBS -l walltime=2400:00:00
    cd $PBS_O_WORKDIR
    #{yield}
  SH
  
  File.open("PBS_#{size}_#{algorithm}.sh", "w") do |file|
    file.write(content.gsub(/^\s*/, ""))
  end
end

def make_pbs_files
  20.step(300, 20).inject({}) do |hash, size| 
    hash.merge(size => (case size; when 1...200 then 100; else 10; end))
  end.each do |size, iterations|
    iterations.times do |i|
      sequence = size.times.inject("") { |string, _| string + %w[A U C G][rand(4)] }
      
      mfe_content = <<-STR
        >
        #{sequence}
        #{ViennaRna::Fold.run(sequence).structure}
      STR
    
      File.open("mfe_size_#{size}_%03d.fa" % (i + 1), "w") do |file|
        file.write(mfe_content.gsub(/^\s*/, ""))
      end
      
      empty_content = <<-STR
        >
        #{sequence}
        #{'.' * size}
      STR
    
      File.open("empty_size_#{size}_%03d.fa" % (i + 1), "w") do |file|
        file.write(empty_content.gsub(/^\s*/, ""))
      end
    end

    # Run each experiment (.sh file) in triplicate
    (1..3).each do |i|
      make_command_file("rnabor_3_mfe", size)   { Dir["mfe_size_#{size}_*.fa"  ].map { |file| "time ./RNAbor -nodangle #{file}" }.join("\n") }
      make_command_file("rnabor_3_empty", size) { Dir["empty_size_#{size}_*.fa"].map { |file| "time ./RNAbor -nodangle #{file}" }.join("\n") }
      
      make_command_file("fftbor_#{i}_mfe", size)   { Dir["mfe_size_#{size}_*.fa"  ].map { |file| "time ./FFTbor -nodangle #{file}" }.join("\n") }
      make_command_file("fftbor_#{i}_empty", size) { Dir["empty_size_#{size}_*.fa"].map { |file| "time ./FFTbor -nodangle #{file}" }.join("\n") }
    end
  end
end