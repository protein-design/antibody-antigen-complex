
CREATE DATABASE COMPLEX;
USE COMPLEX;

CREATE TABLE pdb(
    pdb_id VARCHAR(10),
    pdb_url VARCHAR(50),
    local_path VARCHAR(100),
    release_date DATE,
    resolution float,
    name VARCHAR(500),
    CONSTRAINT pk_pdb PRIMARY KEY (pdb_id)
);

CREATE TABLE compound(
    pdb_id VARCHAR(10),
    compound_no INT,
    chain VARCHAR(20),
    chain_type VARCHAR(20),
    compound_text Text,
    CONSTRAINT pk_compound PRIMARY KEY (pdb_id, compound_no),
    CONSTRAINT fk_compound_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id)
        ON DELETE CASCADE
);


