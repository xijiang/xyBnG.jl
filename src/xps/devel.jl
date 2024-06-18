function devel(; op = 1, nrpt = 1)
    if op == 1
        dosblup(nrpt = nrpt)
    else
        plan = Plan(25, 50, 200)
        println(plan)
        cattle = Cattle(5_000)
        println(cattle)
        milk = Trait("Milk", 0.25, 1000; sex = 0)
        println(milk)
    end
end
