# Run a matlab script
function run_matlab_script(script_path::String, tool_path::String)
    if occursin("matlab", basename(tool_path))
        run(`$tool_path -batch "$script_path"`) # -batch to launch Matlab without graphical interface
    elseif occursin("octave", basename(tool_path))
        run(`$tool_path --no-gui "$script_path"`) # --no-gui to launch Octave without graphical interface
    end
end

# Function to detect if Matlab or Octave is installed
function detect_matlab_or_octave()
    if Sys.iswindows()
        matlab_path = `where matlab`
        octave_path = `where octave`
    else
        matlab_path = `which matlab`
        octave_path = `which octave`
    end
    if success(matlab_path)
        return matlab_path
    elseif success(octave_path)
        return octave_path
    else
        println("Neither Matlab nor Octave where found in your system path. Please provide the path of your file octave-launch.exe or matlab.exe\nIf you want to stop the process and exit, you can write 'exit'")
        println("Example of path: C:/Users/n.barla/AppData/Local/Programs/GNU Octave/Octave-9.2.0/octave-launch.exe")
        println("Before providing tis path, write '1' and press ENTER (this first line will not be read by Julia)")
        tool_path = readline();
        while !is_matlab_octave_path("$tool_path")
            tool_path = readline();
        end
        return tool_path
    end
end

# Function to check if a user provided path correspond to the launcher of Matlab or Octave
function is_matlab_octave_path(tool_path::String)
    if isfile(tool_path) && tool_path[end-3:end] == ".exe"
        if occursin("octave", basename(tool_path)) || occursin("matlab", basename(tool_path))
            return true
        else
            println("Your input file name ($(basename(tool_path))) should be 'octave-launch.exe' or 'matlab.exe'\nPlease retry.")
        end
    elseif occursin("exit", tool_path)
        error("Neither Matlab nor Octave is installed on this system.")
    else
        println("Your input ($tool_path) should be a .exe file or should be 'exit'.\nPlease retry.")
    end
    return false
end
