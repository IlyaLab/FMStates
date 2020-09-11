#!/usr/bin/env ruby
# frozen_string_literal: true

def run(cmdline)
  puts cmdline
  system cmdline
end

def build_images()
  run <<~CMDLINE.strip
    cd docker && \
    docker build --tag fmstates-script:latest --file FMStates.dockerfile .
  CMDLINE
end

def clean_images()
  run <<~CMDLINE.strip
    docker rmi $(docker images -f dangling=true -q)
  CMDLINE
end

def options_users()
  uid = `id -u`.strip
  gid = `id -g`.strip
  "-u #{uid}:#{gid}"
end

def options_logs()
  [
    "max-size=10m",
    "max-file=5",
  ].map { |s| "--log-opt #{s}" }.join(" ")
end

def options_volume()
  "-v `pwd`/project:/project"
end

def options_all()
  <<~CMDLINE.strip
    -p 8888:8888 \
    #{options_logs()} \
    #{options_volume()}
  CMDLINE
end

def run_bash()
  run <<~CMDLINE.strip
    docker run -it --rm --name fmstates-bash \
      #{options_all()} \
      fmstates-script /bin/bash
  CMDLINE
end

def run_notebook()
  run <<~CMDLINE.strip
    docker run -it --rm --name fmstates-notebook \
      #{options_all()} \
      fmstates-script jupyter-notebook
  CMDLINE
end

def stop_all()
  run <<~CMDLINE.strip
    docker stop fmstates-bash fmstates-notebook
  CMDLINE
end

def prepare_directory()
  curr_dir = File.expand_path(".", __dir__)
  Dir.chdir(curr_dir)
end

def main()
  commands = {
    "build" => -> { build_images() },
    "clean" => -> { clean_images() },
    "stop-all" => -> { stop_all() },
    "run-bash" => -> { run_bash() },
    "run-notebook" => -> { run_notebook() },
  }

  if ARGV.empty?
    puts "Available commands:"
    puts commands.keys.join("\n")
    return
  end

  ARGV.each do |name|
    prepare_directory()
    cmd = commands[name]
    if cmd.nil?
      puts "unknown command: #{name}"
    else
      cmd.call()
    end
  end
end

main() if caller.empty?
