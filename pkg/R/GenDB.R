
## we make the driver global!
#drv<-dbDriver("SQLite");


openGenDB <- function(dbName) {
    if(!file.exists(dbName)) cat("Creating a new GenDB.\n")
    else stop("DB already exists")
    con<-dbConnect(drv, dbname = dbName);

    cat("Tables: ")
    print(dbListTables(con))
    con
}

closeGenDB <- function(db) {
    dbDisconnect(db)
}

listGenDB <- function(db) {
    dbListTables(db)
}


# CLOSE THE CONNECTION
#dbDisconnect(con)
#dbUnloadDriver(drv)


