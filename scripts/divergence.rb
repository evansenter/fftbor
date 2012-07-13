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
end

20.step(300, 20).each do |size|
  3.times do
    sequence = size.times.inject("") { |string, _| string + %w[A U C G][rand(4)] }
         
    ["." * sequence.length, ViennaRna::Fold.new(sequence).run.structure].each do |structure|
      rnabor = ViennaRna::Rnabor.run(sequence: sequence, structure: structure)
      fftbor = ViennaRna::Fftbor.run(sequence: sequence, structure: structure)
    
      rnabor_distribution = rnabor.distribution
      fftbor_distribution = fftbor.distribution
    
      # Run.create({
      ap({
        sequence:        sequence, 
        sequence_length: size, 
        structure:       structure, 
        algorithm:       "RNAbor vs. FFTbor (Z_k/Z)", 
        tvd:             Diverge.new(rnabor_distribution, fftbor_distribution).tvd,
        count:           -1,
        fftbor_time:     fftbor.runtime.real,
        rnabor_time:     rnabor.runtime.real
      })
    end
  end
end