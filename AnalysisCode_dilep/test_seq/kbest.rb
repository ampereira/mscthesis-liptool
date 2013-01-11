# !/bin/ruby

require 'date'

class Kbest

# validated
	def initialize(iterations)
		@threshold = 1.05
		@k = 4
		@min_runs = 2
		@max_runs = 4
		@d_iter = iterations
	end

# validated
# Note: I've ensured that the test.sh script runs the
#       app 4 times for the specified number of dilep
#       iterations
	def runapp
		cmd = "./test.sh #{@d_iter}"
		system(cmd)
	end

# Note: the current best value also counts for the
#       k best 
# validated
	def check(best, results)
		mag = best * @threshold
		amount = 0

		results.each do |value|
			if value <= mag
				amount += 1
			end
		end

		if amount < @k
			return false
		else
			return true
		end
	end

# validated
	def print_result(best, flag)
		filename = "kbest_#{@d_iter}.txt"
		date = DateTime.now.to_s

		file = File.open(filename, "a")

		if flag
			string = "#{date}\t#{best}\tvalid\n"
		else
			string = "#{date}\t#{best}\ttimeout\n"
		end
			
		file.write(string)
		file.close
	end

	def run
		i = 1
		best = 10000000000
		results = Array.new
		flag = false

		while i <= @max_runs and !flag do
			runapp
			filename = "time_#{@d_iter}.txt"
			file = File.new(filename,"r")

			while (time = file.gets) do
				a = Float(time)

				if a < best
					best = a
				end

				results << a
			end

# is it in the k best?
			if i >= @min_runs
				flag = check(best, results)
			end

			i += 1
			File.delete(filename)
		end
# always prints the best result even if it has
# no k values within the threshold, in case of
# timeout
		print_result(best, flag)
	end
end

# Run the scritp!
if __FILE__ == $0
	kb = Kbest.new ARGV.first
	kb.run
end

