using SymbolicWedderburn
using PermutationGroups
import SymbolicWedderburn as SW


# Constructing SL(2,3) or binary tetrahedral group as a permutation group of 8 elements
gen1 = perm"(1,2,3,4)(5,6,7,8)" 
gen2 = perm"(1,7,3,5)(2,6,4,8)"
gen3 = perm"(2,6,7)(4,8,5)"  
MyGroup = PermGroup([gen1, gen2, gen3]) # SL(2,3)

# Compute the character table
tbl = SW.CharacterTable(Rational{Int}, MyGroup)

# Get irreducible characters
irreducible_chars = irreducible_characters(tbl)

# Define multiplicities (for simplicity, use twos, it should be even for the quaternion case) 
multiplicities = fill(2, length(irreducible_chars))

# Get real irreducible characters and their multiplicities
real_irreps, real_mults = SW.affordable_real(irreducible_chars, multiplicities)
# Print the Frobinus-Schur indicator of each complex irrep
using SymbolicWedderburn.Characters
for (i, χ) in enumerate(irreducible_characters(tbl))
    fs_indicator = Characters.frobenius_schur(χ)
    println("Frobenius-Schur indicator of $χ: $fs_indicator")
    println()
end
# Print characters of real irreps
println("Real Irreducible Characters:")
for irrep in real_irreps
    println(irrep)
end
# Expected Result: χ does not change if its FS is 1. χ is doubled if its FS is -1. χ is added by another character if its FS is 0.
