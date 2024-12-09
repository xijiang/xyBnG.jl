"""
    savesum(file::AbstractString, df::DataFrame)
Save the summary DataFrame `df` to file `file`. If the file exists, append the
DataFrame to the existing one.
"""
function savesum(file::AbstractString, df::DataFrame)
    if isfile(file)
        tsm = deserialize(file)
        append!(tsm, df)
        serialize(file, tsm)
    else
        serialize(file, df)
    end
end
