function devel()
    @info "DOS OCS"
    lmp = deserialize("rst/tale/founder.lmp")
    irm("rst/tale/ablup.xy", lmp.chip, 1:3200)
end
