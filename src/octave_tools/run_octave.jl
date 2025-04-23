# Run an octave script
function run_octave_script(script_path::String, tool_path::String, kwargs::Dict)
    commands = "run('$script_path');"
    for (key, value) in kwargs
        commands = "$key = '$value'; $commands"
    end
    if occursin("octave", basename(tool_path))
        run(`$tool_path --no-gui --eval "is_octave = '1'; $commands"`) # --no-gui to launch Octave without graphical interface
    else
        @assert occursin("matlab", basename(tool_path)) "The file name of your matlab/octave launcher should contain 'matlab' or octave'\n$tool_path"
        run(`$tool_path -batch  --eval "is_octave = '0'; $commands"`) # FIXME This code has never been tested successfully
    end
end
# Function to detect if Matlab or Octave is installed
function detect_octave()
    if Sys.iswindows()
        get_octave = `where octave`
    else
        get_octave = `which octave`
    end
    if success(get_octave)
        dir = dirname(read(get_octave, String))
        if isfile(joinpath(dir, "octave.exe"))
            return joinpath(dir, "octave.exe")
        elseif isfile(joinpath(dir, "octave-launch.exe"))
            return joinpath(dir, "octave-launch.exe")
        end
    end
    println("Octave was not found in your system path. Please provide the path of your file octave-launch.exe or matlab.exe\nIf you want to stop the process and exit, you can write 'exit'")
    println("To stop having this message, you can add this path to your environment variables")
    println("Example of path: C:/Users/n.barla/AppData/Local/Programs/GNU Octave/Octave-9.2.0/octave-launch.exe")
    println("Before providing this path, write '1' and press ENTER (this first line will not be read by Julia)")
    tool_path = readline();
    while !is_octave_path("$tool_path")
        tool_path = readline();
    end
    return tool_path
end

# Function to check if a user provided path correspond to the launcher of Matlab or Octave
function is_octave_path(tool_path::String)
    if isfile(tool_path) && (tool_path[end-3:end] == ".exe")
        if occursin("octave", basename(tool_path))
            return true
        elseif occursin("matlab", basename(tool_path))
            return true
        else
            println("Your input file name ($(basename(tool_path))) should be 'octave-launch.exe' or 'matlab.exe'\nPlease retry.")
        end
    elseif occursin("exit", tool_path)
        error("Octave is not installed on this system.")
    else
        println("Your input ($tool_path) should be a .exe file or 'exit'.\nPlease retry.")
    end
    return false
end
