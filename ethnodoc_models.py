import os
from sqlalchemy import create_engine, Column, Integer, String, Text, ForeignKey, Table
from sqlalchemy.orm import declarative_base, relationship, sessionmaker

Base = declarative_base()

# Association tables for many-to-many relationships
species_source_association = Table(
    'species_source', Base.metadata,
    Column('species_id', String, ForeignKey('species.species_id')),
    Column('source_id', String, ForeignKey('historical_source.source_id'))
)

claim_compound_association = Table(
    'claim_compound', Base.metadata,
    Column('claim_id', String, ForeignKey('historical_claim.claim_id')),
    Column('compound_id', String, ForeignKey('phytochemical.compound_id')),
    Column('evidence_type', String),
    Column('evidence_strength', String)
)

compound_target_association = Table(
    'compound_target', Base.metadata,
    Column('compound_id', String, ForeignKey('phytochemical.compound_id')),
    Column('target_id', String, ForeignKey('molecular_target.target_id')),
    Column('evidence_type', String)
)

class Species(Base):
    __tablename__ = 'species'
    
    species_id = Column(String, primary_key=True)
    english_name = Column(String)
    scientific_name = Column(String, nullable=False)
    chinese_name = Column(String)
    taxonomy_genus = Column(String)
    
    # Relationships
    historical_sources = relationship("HistoricalSource", secondary=species_source_association, back_populates="species")
    claims = relationship("HistoricalClaim", back_populates="species")
    images = relationship("SpeciesImage", back_populates="species")
    
    def __repr__(self):
        return f"<Species(id='{self.species_id}', scientific='{self.scientific_name}')>"

class SpeciesImage(Base):
    __tablename__ = 'species_image'
    
    image_id = Column(String, primary_key=True)
    species_id = Column(String, ForeignKey('species.species_id'))
    
    filename = Column(String, nullable=False)
    source_name = Column(String)
    source_url = Column(String)
    license_info = Column(String)
    acquisition_date = Column(String)
    
    # Relationships
    species = relationship("Species", back_populates="images")



class HistoricalSource(Base):
    __tablename__ = 'historical_source'
    
    source_id = Column(String, primary_key=True)
    title = Column(String, nullable=False)
    author = Column(String)
    historical_period = Column(String)
    provenance = Column(String)
    
    # Relationships
    species = relationship("Species", secondary=species_source_association, back_populates="historical_sources")
    claims = relationship("HistoricalClaim", back_populates="source")

class HistoricalClaim(Base):
    __tablename__ = 'historical_claim'
    
    claim_id = Column(String, primary_key=True)
    species_id = Column(String, ForeignKey('species.species_id'))
    source_id = Column(String, ForeignKey('historical_source.source_id'))
    
    claim_text = Column(Text)
    translation = Column(Text)
    disease_condition = Column(String)
    
    # Relationships
    species = relationship("Species", back_populates="claims")
    source = relationship("HistoricalSource", back_populates="claims")
    preparations = relationship("Preparation", back_populates="claim")
    compounds = relationship("Phytochemical", secondary=claim_compound_association, back_populates="claims")

class Preparation(Base):
    __tablename__ = 'preparation'
    
    preparation_id = Column(String, primary_key=True)
    claim_id = Column(String, ForeignKey('historical_claim.claim_id'))
    
    plant_part = Column(String)
    processing_method = Column(Text)
    provenance = Column(String)
    
    # Relationships
    claim = relationship("HistoricalClaim", back_populates="preparations")

class Phytochemical(Base):
    __tablename__ = 'phytochemical'
    
    compound_id = Column(String, primary_key=True)
    compound_name = Column(String, nullable=False)
    smiles = Column(String)
    molecular_formula = Column(String)
    molecular_weight = Column(String)
    chemical_class = Column(String)
    provenance = Column(String)
    
    # Relationships
    claims = relationship("HistoricalClaim", secondary=claim_compound_association, back_populates="compounds")
    targets = relationship("MolecularTarget", secondary=compound_target_association, back_populates="compounds")

class MolecularTarget(Base):
    __tablename__ = 'molecular_target'
    
    target_id = Column(String, primary_key=True)
    protein_name = Column(String)
    gene_identifier = Column(String)
    structure_identifier = Column(String) # e.g. PDB ID
    organism = Column(String)
    
    # Relationships
    compounds = relationship("Phytochemical", secondary=compound_target_association, back_populates="targets")
    docking_experiments = relationship("DockingExperiment", back_populates="target")

class DockingExperiment(Base):
    __tablename__ = 'docking_experiment'
    
    experiment_id = Column(String, primary_key=True)
    compound_id = Column(String, ForeignKey('phytochemical.compound_id'))
    target_id = Column(String, ForeignKey('molecular_target.target_id'))
    
    software_version = Column(String)
    docking_score = Column(String)
    poses_json = Column(Text) # Stores the 9 binding modes (Affinity, RMSD l.b., RMSD u.b.)
    random_seed = Column(String)
    parameters_json = Column(Text)
    pose_file_path = Column(String)
    timestamp = Column(String)
    
    # Relationships
    target = relationship("MolecularTarget", back_populates="docking_experiments")

# --- Database Initialization ---
def init_db(db_path="sqlite:///ethnodoc.db"):
    engine = create_engine(db_path, echo=False)
    Base.metadata.create_all(engine)
    return engine

def get_session(engine):
    Session = sessionmaker(bind=engine)
    return Session()
