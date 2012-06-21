require "rubygems"
require "mysql2"
require "active_record"
require "awesome_print"
require "vienna_rna"
require "diverge"

class Object; def this; self; end; end
module Enumerable; def sum; inject(&:+); end; end

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
      table.float   :js_divergence
      table.float   :fftbor_time
      table.float   :rnabor_time
      table.float   :fftbor_partition
      table.float   :rnabor_partition
      table.timestamps
    end 
  end
end

unless ActiveRecord::Base.connection.execute("show tables").map(&:this).flatten.include?("runs")
  BuildRun.up
end

class Run < ActiveRecord::Base
  validates_presence_of :sequence, :sequence_length, :structure, :js_divergence, :fftbor_time, :rnabor_time
  
  def self.generate_sequence(sequence_length)
    sequence_length.times.inject("") { |string, _| string + %w[A U C G][rand(4)] }
  end
end

def normalize(array)
  array.map { |i| i + ((1 - array.sum) / array.length) }
end

ViennaRna.debug = false

20.step(300, 10).each do |size|
  sequence = Run.generate_sequence(size)
         
  ["." * sequence.length, ViennaRna::Fold.new(sequence).run.structure].each do |structure|
    fftbor = ViennaRna::Fftbor.new(sequence: sequence, structure: structure).run
    rnabor = ViennaRna::Rnabor.new(sequence: sequence, structure: structure).run
    
    puts (fftbor_distribution = normalize(fftbor.distribution)).inspect
    puts (rnabor_distribution = normalize(rnabor.distribution)).inspect
    
    attributes = {
      sequence:         sequence, 
      sequence_length:  size, 
      structure:        structure, 
      js_divergence:    Diverge.new(fftbor_distribution, rnabor_distribution).js,
      fftbor_time:      fftbor.runtime.real,
      rnabor_time:      rnabor.runtime.real,
      fftbor_partition: fftbor.partition,
      rnabor_partition: rnabor.partition,
    }
    
    ARGV.last == "save" ? Run.create(attributes) : ap(attributes)
  end
end