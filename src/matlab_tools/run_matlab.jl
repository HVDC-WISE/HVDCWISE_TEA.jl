# Run a matlab script
function run_matlab_script(script_path::String, tool_path::String)
    if occursin("matlab", basename(tool_path))
        run(`$tool_path -batch "$script_path"`) # -batch to launch Matlab without graphical interface
    elseif occursin("octave", basename(tool_path))
        run(`$tool_path --no-gui "$script_path"`) # --no-gui to launch Octave without graphical interface
    end
end
# run(`"C:\\Users\\n.barla\\AppData\\Local\\Programs\\GNU Octave\\Octave-9.2.0\\octave.vbs\r\n" --no-gui 'C:\\Users\\n.barla\\Documents\\Local_codes\\HVDCWISE_TEA.jl\\src\\matlab_tools\\build_availability_series.m'`)
# Function to detect if Matlab or Octave is installed
function detect_matlab_or_octave()
    if Sys.iswindows()
        get_matlab = `where matlab`
        get_octave = `where octave`
    else
        get_matlab = `which matlab`
        get_octave = `which octave`
    end
    if success(get_matlab)
        dir = dirname(read(get_matlab, String))
        if isfile(joinpath(dir, "matlab.exe"))
            return joinpath(dir, "matlab.exe")
        end
    elseif success(get_octave)
        dir = dirname(read(get_octave, String))
        if isfile(joinpath(dir, "octave.exe"))
            return joinpath(dir, "octave.exe")
        elseif isfile(joinpath(dir, "octave-launch.exe"))
            return joinpath(dir, "octave-launch.exe")
        end
    end
    println("Neither Matlab nor Octave where found in your system path. Please provide the path of your file octave-launch.exe or matlab.exe\nIf you want to stop the process and exit, you can write 'exit'")
    println("To stop having this message, you can add this path to your environment variables")
    println("Example of path: C:/Users/n.barla/AppData/Local/Programs/GNU Octave/Octave-9.2.0/octave-launch.exe")
    println("Before providing this path, write '1' and press ENTER (this first line will not be read by Julia)")
    tool_path = readline();
    while !is_matlab_octave_path("$tool_path")
        tool_path = readline();
    end
    return tool_path
end

# Function to check if a user provided path correspond to the launcher of Matlab or Octave
function is_matlab_octave_path(tool_path::String)
    if isfile(tool_path) && (tool_path[end-3:end] == ".exe" || tool_path[end-3:end] == ".vbs")
        if occursin("octave", basename(tool_path)) || occursin("matlab", basename(tool_path))
            return true
        else
            println("Your input file name ($(basename(tool_path))) should be 'octave-launch.exe' or 'matlab.exe'\nPlease retry.")
        end
    elseif occursin("exit", tool_path)
        error("Neither Matlab nor Octave is installed on this system.")
    else
        println("Your input ($tool_path) should be a .exe file, a .vbs file, or 'exit'.\nPlease retry.")
    end
    return false
end
