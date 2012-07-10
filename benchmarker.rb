require "rubygems"
require "vienna_rna"
require "awesome_print"
require "resque"
require "./benchmark_run.rb"
require "./benchmark_job.rb"

class Object; def this; self; end; end

Run.connect

class BuildRun < ActiveRecord::Migration
  def self.up
    create_table :runs do |table|
      table.string  :sequence
      table.string  :structure
      table.integer :sequence_length
      table.string  :algorithm
      table.decimal :time, precision: 20, scale: 3
      table.timestamps
    end 
  end
end

unless ActiveRecord::Base.connection.execute("show tables").map(&:this).flatten.include?("runs")
  BuildRun.up
end

ViennaRna.debug = false

# [20, 40, 60, 80, 100, 120, 140, 160, 200, 250, 300].inject({}) do |hash, size| 
#   hash.merge(size => (case size; when 1..160 then 100; else 3; end))
# end.each do |size, iterations|
#   iterations.times do |i|
#     sequence = Run.generate_sequence(size)
#     
#     Resque.enqueue(BenchmarkJob, { algorithm: :rnabor, sequence: sequence })
#     Resque.enqueue(BenchmarkJob, { algorithm: :fftbor, sequence: sequence })
#   end
# end

# Run.find_in_batches do |runs|
#   runs.each do |run|
#     Resque.enqueue(DuplicateJob, { id: run.id })
#     Resque.enqueue(DuplicateJob, { id: run.id })
#   end
# end