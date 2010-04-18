#ifndef GENERATOR_HARMONIC_H
#define GENERATOR_HARMONIC_H

Generator * harmonic_generator_creator();
static int harmonic_generator_creator_registred = register_generator("harmonic", harmonic_generator_creator);

#endif /* GENERATOR_HARMONIC_H */
