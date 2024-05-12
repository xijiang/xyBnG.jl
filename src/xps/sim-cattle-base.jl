"""
    function sim_cattle_base()

This sample function is to simulate a cattle base population, using `stdpopsim`.
"""
function sim_cattle_base()
    cattle = Cattle("cattle", 5_000)
    sim_base(cattle, "rst/base")
end
