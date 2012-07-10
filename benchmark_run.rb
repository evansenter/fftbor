require "mysql2"
require "active_record"

class Run < ActiveRecord::Base
  validates_presence_of :sequence, :sequence_length, :algorithm
  
  def self.connect
    ActiveRecord::Base.establish_connection(config = { adapter: "mysql2", username: "root", reconnect: true })

    unless ActiveRecord::Base.connection.execute("show databases").map { |i| i }.flatten.include?("fftbor_performance")
      ActiveRecord::Base.connection.create_database("fftbor_performance")
    end

    ActiveRecord::Base.establish_connection(config.merge(database: "fftbor_performance"))
    
    inline_rails if defined?(inline_rails)
  end
  
  def self.generate_sequence(sequence_length)
    sequence_length.times.inject("") { |string, _| string + %w[A U C G][rand(4)] }
  end
end