require "resque"
require "awesome_print"
require "vienna_rna"
require "./benchmark_run.rb"

module BenchmarkJob
  @queue = :benchmarking

  def self.perform(params)
    Run.connect
    
    size = params["sequence"].length
    
    ["." * size, ViennaRna::Fold.run(params["sequence"]).structure].each do |structure|
      results = case params["algorithm"]
      when "rnabor" then ViennaRna::Rnabor.new(sequence: params["sequence"], structure: structure).run
      when "fftbor" then ViennaRna::Fftbor.new(sequence: params["sequence"], structure: structure).run
      end
      
      Run.create({
        sequence:        params["sequence"], 
        structure:       structure, 
        sequence_length: size, 
        algorithm:       params["algorithm"], 
        time:            results.runtime.real
      })
    end
  end
end

module DuplicateJob
  @queue = :benchmarking

  def self.perform(params)
    Run.connect
    
    run = Run.find(params["id"])
    
    results = case run.algorithm
    when "rnabor" then ViennaRna::Rnabor.new(sequence: run.sequence, structure: run.structure).run
    when "fftbor" then ViennaRna::Fftbor.new(sequence: run.sequence, structure: run.structure).run
    end
      
    Run.create({
      sequence:        run.sequence, 
      structure:       run.structure, 
      sequence_length: run.sequence.length, 
      algorithm:       run.algorithm, 
      time:            results.runtime.real
    })
  end
end