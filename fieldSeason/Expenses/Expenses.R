budget <- read.csv("/Users/laurenponisio/Documents/Sky Islands/Expenses/Expenses.csv")

aggregate(budget$amount, list(budget$expense), sum)
total <- sum(budget$amount)

food <- sum(budget$amount[budget$expense == "food"])/34

notfood <- sum(budget$amount[!budget$expense == "food"])