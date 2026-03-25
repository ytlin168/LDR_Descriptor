use strict;
use warnings;
use MdmCommands;
use MdmDiscoveryScript;
use SdmDiscoveryScript;
use ForceFieldDiscoveryScript;
use ChartDiscoveryScript;
use Text::CSV;
use File::Basename;
use File::Path qw(make_path);

# Define document name and target protein name
my $doc_name = "3B1M_Diversity_Flexibility";  # Name of the molecular window to analyze
my $protein_name = "3B1M";  # Name of the target protein
my $Hbond_length = "10";  # Hydrogen bond length
my $test_name = $doc_name . "(" .$Hbond_length. ")";
my $save_path = "C:/Users/ytlin/Desktop/results/receptor/GOLD_diversity_0.5/0.9_0.45/" . $test_name . "/";  # Save path for each ligand configuration (target folder must be created manually)
my $property = 'Gold.Goldscore.Fitness';  # Name of the scoring function

# Set the save path for the hydrogen bond configuration summary results
my $output_folder = 'C:/Users/ytlin/Desktop/perldealed/';

# Read all documents (all currently open windows in DS)
my $Documents = DiscoveryScript::Documents();
my $doc;

# Iterate through windows to find the target window
for (my $i = 0; $i <= $Documents->Count; $i++) {
    my $document;
    eval { $document = $Documents->Item($i); };
    if ($@) {
        print "Error accessing document item: $@\n";
        next;
    }

    # Attempt to read the window name
    my $Name;
    eval { $Name = $document->ModelName; };
    if ($@) {
        print "Error retrieving document ModelName: $@\n";
    } else {
        if (defined $Name && $Name eq $doc_name) {       
            $doc = $document;
            last;
        }
    }
}

# Verify if the target document was successfully loaded
if (defined $doc) {
    print "Document found.\n";

    # Read target protein
    my $Molecules = $doc->AllMolecules;
    if (!$Molecules) {
        print "Error: No molecules found in document.\n";
        exit;
    }

    my $protein;
    for (my $i = 0; $i < $Molecules->Count; $i++) {
        my $mol;
        eval { $mol = $Molecules->Item($i); };
        if ($@) {
            print "Error accessing molecule item: $@\n";
            next;
        }

        if (defined $mol->Name && $mol->Name eq $protein_name) {
            $protein = $mol;
            last;
        }
    }

    # Verify if the target protein was successfully loaded
    if (defined $protein) {
        print "Protein $protein_name found.\n";

        # Read target amino acids (residues)
        my $Residues = $protein->Residues;
        my %target_residues = (
            "HIS323" => "HIS323323",
            "HIS449" => "HIS449449",
            "TYR473" => "TYR473473",
        );

        my @residues_to_move;

        if ($Residues) {
            for (my $i = 0; $i < $Residues->Count; $i++) {
                my $residue;
                eval { $residue = $Residues->Item($i); };
                if ($@) {
                    print "Error accessing residue item: $@\n";
                    next;
                }
          
                my $residue_name = $residue->Name . $residue->Id;
                if (exists $target_residues{$residue->Name}) {
                    push @residues_to_move, { residue => $residue, new_name => $target_residues{$residue->Name} };
                    $target_residues{$residue->Name} = 1;
                }
            }
        }

        # Verify if all target residues were found
        my $all_residues_found = 1;
        for my $key (keys %target_residues) {
            if ($target_residues{$key} ne "1") {
                $all_residues_found = 0;
                print "Error: Residue $key was not found.\n";
            }
        }
        # Isolate the target residues from the protein
        if ($all_residues_found) {
            foreach my $res_info (@residues_to_move) {
                eval {
                    my $new_molecule = CopyIntoNewMolecule($doc, $res_info->{residue});
                    $new_molecule->Name($res_info->{new_name});
                    for (my $i = 0; $i < $new_molecule->Chains->Count; $i++){
                        eval {$new_molecule->Chains->Item($i)->Name($res_info->{new_name});};}
                };
                if ($@) {
                    print "Error moving residue to new molecule: $@\n";
                    next;
                }
                print "Residue moved to new molecule and renamed to ", $res_info->{new_name}, ".\n";
            }
            
            if (defined $doc) {
                print "Document found.\n";
                print "There are " . $Molecules->Count . " molecules in the file.\n";
            }
           
            # Delete the protein that no longer contains the target residues
            eval { Delete($doc, $protein); };
            if ($@) {
                print "Error deleting original protein: $@\n";
            } else {
                print "Original protein deleted.\n";
            }

            # Save the modified molecular window
            my $dealed_path = $save_path . "dealed1.dsv";  # Set to desired filename using the previously defined path
            $doc->Save($dealed_path, 'dsv');
            print "The modified document has been saved.\n";

        } else {
            print "Error: Not all specified residues were found.\n";
        }

        # Isolate each ligand's docking pose into a new window; perform the following functions for each isolated ligand
        my %ligand_groups;
        for (my $i = 0; $i < $Molecules->Count; $i++) {
            my $ligand;
            eval { $ligand = $Molecules->Item($i); };
            if ($@) {
                print "Error accessing molecule item: $@\n";
                next;
            }
            next if $ligand->Name eq $protein_name;

            push @{$ligand_groups{$ligand->Name}}, $ligand;
        }

        foreach my $ligand_name (keys %ligand_groups) {
            my $new_doc = Mdm::Document::Create();
            foreach my $ligand (@{$ligand_groups{$ligand_name}}) {
                eval {
                    my $ligand_copy = $ligand->Copy();
                    my $newDocumentRoot = $new_doc->RootObject;
                    my $new_molecule = $new_doc->ChangeObjectParent($ligand_copy, $newDocumentRoot);
                    Delete($doc, $ligand);
                };     
                if ($@) {
                    print "Error copying ligand to new document: $@\n";
                    next;
                }
            }

            # Ensure each new window contains the three critical residues: HIS323, HIS449, and TYR473
            my %target_residue_names = map { $_ => 1 } ("HIS323323", "HIS449449", "TYR473473");
            foreach my $ligand_residue_name ("HIS323323", "HIS449449", "TYR473473") {
                my $ligand_residue;
                my $Molecules2 = $doc->AllMolecules;
                eval {
                    for (my $i = 0; $i < $Molecules2->Count; $i++) {
                        my $mol;
                        eval { $mol = $Molecules2->Item($i); };
                        if ($@) {
                            print "Error accessing molecule item: $@\n";
                            next;
                        }
                        if ($mol->Name eq $ligand_residue_name) {
                            $ligand_residue = $mol;                            
                            last;
                        }
                    }                
                };
                if ($@) {
                    print "Error accessing molecule $ligand_residue_name: $@\n";
                    next;
                }
                if (defined $ligand_residue) {
                    eval {
                       my $ligand_residue_copy = $ligand_residue->Copy();
                       my $newDocumentRoot = $new_doc->RootObject;
                       my $new_molecule = $new_doc->ChangeObjectParent($ligand_residue_copy, $newDocumentRoot);  
                    };
                    if ($@) {
                        print "Error copying ligand $ligand_residue_name to new document: $@\n";
                        next;
                    }
                } else {
                    print "Ligand $ligand_residue_name not found.\n";
                }
            }

            # Rename other configurations based on Index sequence, excluding HIS323, HIS449, and TYR473
            my $index = 1;
            for (my $i = 0; $i < $new_doc->AllMolecules->Count; $i++) {
                my $ligand = $new_doc->AllMolecules->Item($i);
                $ligand->SetVisible(True);
                unless (exists $target_residue_names{$ligand->Name}) {
                    $ligand->Name($index);
                    $index++;
                }
            }

            # Read the $property attribute (defined scoring function) for each pose and sort them
            my @ligands_with_score;
            for (my $i = 0; $i < $new_doc->AllMolecules->Count; $i++) {
                my $ligand = $new_doc->AllMolecules->Item($i);
                my $score = eval { $ligand->GetProperty($property) };
                if ($@) {
                    print "Error accessing property: $@\n";
                    $score = undef;  # If score retrieval fails, set score to undef
                }
                push @ligands_with_score, { ligand => $ligand, score => $score } if defined $score && $score ne '';
            }
            @ligands_with_score = sort { $b->{score} <=> $a->{score} } @ligands_with_score;

            # Select only the top 20 poses plus HIS323, HIS449, and TYR473
            my @selected_ligands;
            %target_residue_names = map { $_ => 1 } ("HIS323323", "HIS449449", "TYR473473");
            for (my $i = 0; $i < 20 && $i < @ligands_with_score; $i++) {
                push @selected_ligands, $ligands_with_score[$i];
            }
            foreach my $ligand_info (@ligands_with_score) {
                if (exists $target_residue_names{$ligand_info->{ligand}->Name}) {
                    push @selected_ligands, $ligand_info unless grep { $_->{ligand}->Name eq $ligand_info->{ligand}->Name } @selected_ligands;
                }
            }

            # List the top 20 names and scores in the output window
            print "Top 20 ligands and their scores:\n";
            for (my $i = 0; $i < @selected_ligands; $i++) {
                print $selected_ligands[$i]->{ligand}->Name . ": " . $selected_ligands[$i]->{score} . "\n";
            }

              # Keep selected ligands and delete the others
            my %selected_ligands_hash = map { $_->{ligand}->Name => 1 } grep { defined $_->{ligand} } @selected_ligands;
            my $molecules_to_keep = $new_doc->AllMolecules;
            for (my $i = $molecules_to_keep->Count - 1; $i >= 0; $i--) {
                my $ligand = $molecules_to_keep->Item($i);
                unless (exists $selected_ligands_hash{$ligand->Name} || exists $target_residue_names{$ligand->Name}) {
                    eval { Delete($new_doc, $ligand); };
                    if ($@) {
                        print "Error deleting ligand: $@\n";
                    } else {
                        print "Deleted ligand: ", $ligand->Name, "\n";
                    }
                }
            }
            
            # Create Hbond Monitor to predict hydrogen bonds
            my $residue_atom;
            my $ligand_residue;
            foreach my $ligand_residue_name ("HIS323323", "HIS449449", "TYR473473") {                
                my $Molecules2 = $new_doc->Molecules;
                eval {
                    for (my $i = 0; $i < $Molecules2->Count; $i++) {
                        my $mol;
                        eval { $mol = $Molecules2->Item($i); };
                        if ($@) {
                            print "Error accessing molecule item: $@\n";
                            next;
                        }
                        if ($mol->Name eq $ligand_residue_name) {
                            $ligand_residue = $mol;
                            $ligand_residue->Select();
                            last;
                        }else {
                           print "Error: Ligand residue $ligand_residue_name not found.\n";}
                      }
                    }
                };
                                
            
            
            if (defined $ligand_residue) {
                    $residue_atom = $new_doc->SelectedObjects->Atoms;}
            # Configure Hbond Monitor settings        
            if ($residue_atom) {
                my $Hbondmonitor = $new_doc->CreateNonbondMonitor($residue_atom, Mdm::interactionsExternalToMonitoredScope, Mdm::intermolecularInteractions);
                $Hbondmonitor->Name('ligand_residue_Hmonitor');
                $Hbondmonitor->NonbondCriterion(Mdm::hydrogenBond_Strong_DonorAcceptor_DistanceMaximum) = $Hbond_length;  # Strong H-bond length
                $Hbondmonitor->NonbondCriterion(Mdm::hydrogenBond_Weak_DonorAcceptor_DistanceMaximum) = $Hbond_length;   # Weak H-bond length
                $new_doc->UpdateObjects;
            
               # Export prediction results to a .CSV file
                my $non_bonds = $Hbondmonitor->Nonbonds;
                my $csv_path = $save_path . $ligand_name . "_non_bond_data.csv";
                open(my $fh, '>', $csv_path) or die "Could not open file '$csv_path' $!";
                print $fh "Name,ID,Visible,Visibility Locked,Color,Parent,XYZ:X,XYZ:Y,XYZ:Z,Distance,Category,Types,From,From Chemistry,To,To Chemistry,Angle DHA,Angle HAY,Angle XDA,Angle DAY\n";  # Add Excel headers based on attribute names
                
                for (my $i = 0; $i < $non_bonds->Count; $i++) {  # Starting index modified to 0; adjusted loop condition
                    eval {
                        my $interaction = $non_bonds->Item($i);

                    my $name = eval { $interaction->Name } // 'N/A';
                    my $id = eval { $interaction->ID } // 'N/A';
                    my $visible = eval { $interaction->GetProperty("Visible") } // 'N/A';
                    my $visibility_locked = eval { $interaction->GetProperty("Visibility Locked") } // 'N/A';
                    my $color = eval { $interaction->GetProperty("Color") } // 'N/A';
                    my $parent = eval { $interaction->GetProperty("Parent") } // 'N/A';
                    my $xyz_x = eval { $interaction->GetProperty("XYZ:X") } // 'N/A';
                    my $xyz_y = eval { $interaction->GetProperty("XYZ:Y") } // 'N/A';
                    my $xyz_z = eval { $interaction->GetProperty("XYZ:Z") } // 'N/A';
                    my $distance = eval { $interaction->GetProperty("Distance") } // 'N/A';
                    my $category = eval { $interaction->GetProperty("Category") } // 'N/A';
                    my $types = eval { $interaction->GetProperty("Types") } // 'N/A';
                    my $from = eval { $interaction->GetProperty("From") } // 'N/A';
                    my $from_chemistry = eval { $interaction->GetProperty("From Chemistry") } // 'N/A';
                    my $to = eval { $interaction->GetProperty("To") } // 'N/A';
                    my $to_chemistry = eval { $interaction->GetProperty("To Chemistry") } // 'N/A';
                    my $angle_dha = eval { $interaction->GetProperty("Angle DHA") } // 'N/A';
                    my $angle_hay = eval { $interaction->GetProperty("Angle HAY") } // 'N/A';
                    my $angle_xda = eval { $interaction->GetProperty("Angle XDA") } // 'N/A';
                    my $angle_day = eval { $interaction->GetProperty("Angle DAY") } // 'N/A';
            
                    print $fh join(",", $name, $id, $visible, $visibility_locked, $color, $parent, $xyz_x, $xyz_y, $xyz_z, $distance, $category, $types, $from, $from_chemistry, $to, $to_chemistry, $angle_dha, $angle_hay, $angle_xda, $angle_day) . "\n";
                };
                if ($@) {
                    print "Error processing interaction at index $i: $@\n";
                }
            }
            
            close $fh;
            print "Non-bond monitor results saved to $csv_path\n";
        }
                                             
                    
            # Save the molecular window after H-bond prediction as a .dsv file
            my $new_file_name = $save_path . $ligand_name . ".dsv";  # Save based on the initial configuration
            $new_doc->Save($new_file_name, 'dsv');
            print "New document for ligand $ligand_name has been saved as $new_file_name.\n";

            $new_doc->Close();
        }

    } else {
        print "Protein $protein_name not found.\n";
    }

    # Close the molecular window
    $doc->Close();
} else {
    print "Document not found.\n";
}

# Hydrogen bond pose organization
# Set the folder path for the .CSV files to be processed
opendir(my $dh, $save_path) or die "Cannot read folder '$save_path': $!";

# Initialize array to store all results
my @all_pose_counts;
my @required_columns = ("Name", "Types", "From", "To");

# Get all .CSV files in the folder
my @file_list = grep { /\.csv$/ && -f "$save_path/$_" } readdir($dh);
closedir $dh;

# Check if the file list is empty
die "No CSV files found in the specified folder." unless @file_list;

# Initialize CSV parser
my $csv = Text::CSV->new({ binary => 1, auto_diag => 1 });

# Read and process each .CSV file
foreach my $file (@file_list) {
    my $file_path = "$save_path/$file";
    open my $fh, "<:encoding(utf8)", $file_path or warn "Cannot open file '$file_path': $!\n" and next;

    # Read headers and verify required columns
    my $headers = $csv->getline($fh);
    my %col_idx;
    @col_idx{@$headers} = (0..$#$headers);

    # Check if mandatory columns exist
    my $missing_cols = grep { !exists $col_idx{$_} } @required_columns;
    if ($missing_cols) {
        warn "File '$file' is missing required columns, skipping processing.";
        close $fh;
        next;
    }

    # Initialize storage for matching residues
    my %unique_residues;
    while (my $row = $csv->getline($fh)) {
        next unless $row->[$col_idx{"Types"}]  =~ /Conventional Hydrogen Bond/;
        foreach my $col ("From", "To") {
            my ($residue) = $row->[$col_idx{$col}] =~ /^([^:]+)/;
            next if $residue =~ /^(HIS449449|TYR473473|HIS323323)$/;
            $unique_residues{$residue}++;
        }
    }
    close $fh;

    # Count unique residues and store the filename
    my $unique_count = scalar keys %unique_residues;
    my $file_name = fileparse($file, qr/_non.*/);
    push @all_pose_counts, [$file_name, $unique_count];
}
# Set .CSV output filename
my $output_file = $output_folder . $test_name . ".csv";

# Create output folder if it doesn't exist
make_path($output_folder) unless -d $output_folder;

# Export results to a new .CSV file
open my $out_fh, ">:encoding(utf8)", $output_file or die "Cannot open file '$output_file': $!";
$csv->say($out_fh, ["File", "Unique_Residue_Count"]);
$csv->say($out_fh, $_) for @all_pose_counts;
close $out_fh;

print "Results saved to $output_file\n";