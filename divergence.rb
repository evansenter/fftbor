require "rubygems"
require "mysql2"
require "active_record"
require "awesome_print"
require "vienna_rna"
require "diverge"

class Object; def this; self; end; end

ActiveRecord::Base.establish_connection(config = { adapter: "mysql2", username: "root", reconnect: true })

unless ActiveRecord::Base.connection.execute("show databases").map(&:this).flatten.include?("fftbor_divergence")
  ActiveRecord::Base.connection.create_database("fftbor_divergence")
end

ActiveRecord::Base.establish_connection(config.merge(database: "fftbor_divergence"))

class BuildRun < ActiveRecord::Migration
  def self.up
    create_table :runs do |table|
      table.string  :sequence
      table.integer :sequence_length
      table.string  :structure
      table.string  :algorithm
      table.float   :tvd
      table.float   :count
      table.float   :fftbor_time
      table.float   :rnabor_time
      table.timestamps
    end 
  end
end

unless ActiveRecord::Base.connection.execute("show tables").map(&:this).flatten.include?("runs")
  BuildRun.up
end

class Run < ActiveRecord::Base
  validates_presence_of :sequence, :sequence_length, :structure, :algorithm, :tvd, :count, :fftbor_time, :rnabor_time
  
  def self.generate_sequence(sequence_length)
    sequence_length.times.inject("") { |string, _| string + %w[A U C G][rand(4)] }
  end
end

20.step(300, 20).each do |size|
  sequence = Run.generate_sequence(size)
         
  ["." * sequence.length, ViennaRna::Fold.new(sequence).run.structure].each do |structure|
    fftbor = ViennaRna::Fftbor.new(sequence: sequence, structure: structure).run
    rnabor = ViennaRna::Rnabor.new(sequence: sequence, structure: structure).run
    
    fftbor_distribution = fftbor.distribution
    rnabor_distribution = rnabor.distribution
    
    Run.create({
      sequence:        sequence, 
      sequence_length: size, 
      structure:       structure, 
      algorithm:       "Structures per shell normalized", 
      tvd:             Diverge.new(fftbor_distribution, rnabor_distribution).tvd,
      count:           rnabor.total_count,
      fftbor_time:     fftbor.runtime.real,
      rnabor_time:     rnabor.runtime.real
    })
  end
end